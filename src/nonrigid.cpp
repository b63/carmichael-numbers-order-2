#include <iostream>
#include <vector>
#include <functional>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <nonrigid.h>
#include <counting_factors.h>
#include <util.h>
#include <config.h>


void
generate_nonrigid_cprimes(long p0, long p1,
        const Factorization &L_fact, const Factorization &M_fact,
        long sieve_size)
{
    NTL::ZZ L, M;
    multiply_factors(L, L_fact.primes, L_fact.powers);
    multiply_factors(M, M_fact.primes, M_fact.powers);

    std::vector<long> primes;

    size_t ind = 0;
    const size_t Lprimes_size = L_fact.primes.size();
    const std::vector<long> &Lprimes = L_fact.primes;
    std::function<bool(size_t, long)> filter = [&L, &ind, &Lprimes, &Lprimes_size] (size_t i, long n) -> bool {
            /* assume that L_fact.primes is sorted, and check 
             * if n is in it */
            if(ind < Lprimes_size)
            {
                for(; n < Lprimes[ind]; ++ind)
                    ; /* NOP */
                if ( ind >= Lprimes_size  || n == Lprimes[ind])
                    return false;
            }
            NTL::ZZ p2 { n };
            NTL::sqr(p2, p2);

            if (NTL::divide(L, p2-1))
                return true;
            return false;
        };

    sieve_primes(primes, sieve_size, &filter);
#if LOG_LEVEL >= 1
    std::cout << "sieze size " << sieve_size << ", " << primes.size() << " primes\n";
#endif

    std::vector<std::vector<size_t> > cprimes;
    subset_product_brute_force(cprimes, primes);
}

void
subset_product_brute_force(std::vector<std::vector<size_t>> &cprimes, const std::vector<long> &primes)
{
    size_t num_primes = primes.size();
    // go through all possible subset sizes starting from 2
    size_t count = 0;
    size_t itrcount = 0;
    for (size_t t = 2; t <= num_primes; t++) 
    {
        size_t factors = t;
        std::vector<size_t> index_stack = {0};
#if LOG_LEVEL >= 1
        std::cout << "checking subsets of size " << t << " ...\n";
#endif

        /* cache of products */
        std::vector<NTL::ZZ> products; 
        products.reserve(t);

        // go through subsets of size factors suing a 'stack' (or odometer?)
        size_t index_size;
        while ((index_size = index_stack.size()) > 0) 
        {
#if LOG_LEVEL >= 2
            printVec<NTL::ZZ>(products);
            printVec<size_t>(index_stack);
            std::cout << "\n";
#endif
            itrcount++;
            const size_t prod_size  = products.size();
            size_t top = index_size - 1;
            size_t prod_top = prod_size - 1;

            if (index_stack[top] < num_primes) 
            {
                if (prod_size+1 < factors) 
                {
                    NTL::ZZ prod { 1 };
                    if (prod_size)
                        prod = products[prod_top];
                    products.push_back( prod * primes[index_stack[top]] );
                }

                if (index_size < factors) 
                {
                    index_stack.push_back(index_stack[top]+1);
                    continue;
                }
                /* IMPORTANT: must have at least TWO factors */
                /*            otherwise segfault  */
                NTL::ZZ prod = products[top - 1] * primes[index_stack[top]];

# if LOG_LEVEL >= 1
                printVec<NTL::ZZ>(products);
                printVec<size_t>(index_stack);
                std::cout << " = " << prod << "\n";
#endif

                count++;
                index_stack[top]++;
            }
            else
            {
                index_stack.pop_back();
                if (top > 0) {
                    const size_t i = ++index_stack[top - 1];
                    products.pop_back();
                    if (i < num_primes) 
                    {
                        NTL::ZZ prod { 1 };
                        if (prod_size > 1)
                            prod = products[prod_top-1];
                        products.push_back(prod * NTL::ZZ{primes[i]});
                    }
                }

                continue;
            }
        }
    }

    std::cout << "total count: " << count << "\n";
    std::cout << "op count: " << itrcount << "\n";
}

void
get_nonrigid_factors(std::vector<std::array<long, 2> > nonrigid_factors, 
        const NTL::ZZ &L_val, const NTL::ZZ &M_val,
        size_t sieve_size)
{
#if LOG_LEVEL >= 1
    std::cout << "sieve size " << sieve_size << ", filtering ...\n";
#endif
    std::vector<long> primes;
    std::function<bool(size_t, long)> filter { [&L_val, &M_val](size_t index, long prime) -> bool {
#if LOG_LEVEL >= 2
            if ((index & STEP_MASK) == 0) std::cerr << std::setw(10) << index << "\r";
#endif
            if (NTL::divide(M_val, prime+1) && NTL::divide(L_val, prime-1))
                return true;
            else
                return false;
        }
    };
    sieve_primes(primes, sieve_size, &filter);
#if LOG_LEVEL >= 2
    std::cerr << std::setw(10) << sieve_size << "\n";
#endif

    const size_t num_primes { primes.size() };
    printVec<long>(primes);
    std::cout << "\n";

#if LOG_LEVEL >= 1
    std::cout << "filtering " << num_primes << " primes for possible nonrigid factors ...\n";
#endif

    for(size_t i { 0 }; i < num_primes; ++i)
    {
        const long p0 { primes[i] };
        const NTL::ZZ p0_zz { p0 };
        const NTL::ZZ p2_1 {NTL::sqr(p0_zz)-1};

        const NTL::ZZ g1 { NTL::GCD(L_val, p2_1) };
        const NTL::ZZ g2 { NTL::GCD(M_val, p2_1) };

        for(size_t j {i+1}; j < num_primes; ++j)
        {
            const long p1 { primes[j] };
            const NTL::ZZ p0p1 { p0_zz * p1 };
            if (NTL::divide(p0p1-1, g1) && NTL::divide(p0p1+1, g2))
                nonrigid_factors.push_back(std::move(std::array<long, 2>{p0, p1}));
        }
#if LOG_LEVEL >= 2
        if((i&STEP_MASK) == 0) std::cerr << std::setw(10) << i << "/" << num_primes << "\r";
#endif
    }

#if LOG_LEVEL >= 2
    std::cerr << std::setw(10) << num_primes << "/" << num_primes << "\n";
#endif
}
