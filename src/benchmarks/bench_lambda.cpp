#include <vector>
#include <functional>
#include <cmath>

#include <NTL/ZZ.h>
#include <benchmark/benchmark.h>
#include <util.h>

void
subset_product_brute_force_no_lambda(const NTL::ZZ &L, const NTL::ZZ &M,
        const long p0, const long p1, const std::vector<long> &primes)
{
    NTL::ZZ z_p0 {p0};
    NTL::ZZ z_p1 {p1};
    const NTL::ZZ p0p1 { p0 * p1 };
    const NTL::ZZ p02_1 {NTL::sqr(p0) -1};
    const NTL::ZZ p12_1 {NTL::sqr(p1) -1};

    size_t num_primes = primes.size();
    // go through all possible subset sizes starting from 2
    size_t count = 0;
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
                NTL::ZZ n0 { products[top - 1] * primes[index_stack[top]] };

                const NTL::ZZ n0_1 {n0-1};
                if (NTL::divide(n0_1, p02_1) && NTL::divide(n0_1, p12_1))
                {
                    NTL::ZZ n {n0 * p0p1};
                    if (NTL::divide(n-1, L) && NTL::divide(n+1, M))
                    {
                        std::cout << n << "\n";
                    }
                }

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

# if LOG_LEVEL >= 1
    std::cout << "total count: " << count << "\n";
#endif
}


void
subset_product_brute_force_lambda(const std::function<bool(NTL::ZZ&)> filter, const std::vector<long> &primes)
{
    size_t num_primes = primes.size();
    // go through all possible subset sizes starting from 2
    size_t count = 0;
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
                NTL::ZZ n0 = products[top - 1] * primes[index_stack[top]];
                benchmark::DoNotOptimize(n0);
                if (filter(n0))
                    std::cout << n0 << "\n";

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

# if LOG_LEVEL >= 1
    std::cout << "total count: " << count << "\n";
#endif
}
