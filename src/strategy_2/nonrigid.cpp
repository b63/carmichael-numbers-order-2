#include <iostream>
#include <vector>
#include <functional>
#include <climits>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <strategy_2/nonrigid.h>
#include <counting_factors.h>
#include <util.h>
#include <config.h>


void
generate_possible_factors(std::vector<std::array<long,2> > &factors, const NTL::ZZ &L_val, const long max)
{
    NTL::PrimeSeq s;
    std::vector<long> primes;

    long p { s.next() };
    while (p != 0 && (!max || p < max))
    {
        primes.push_back(p);
        p = s.next();
    }

    const size_t num_primes { primes.size() };
    for(size_t i { 0 }; i < num_primes; ++i)
    {
        const long p0 { primes[i] };
        const NTL::ZZ p0_zz { primes[i] };
        if (NTL::divide(L_val, p0_zz)) 
            continue;

        const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
        if (NTL::divide(L_val, p02_1)) 
            continue;

        const NTL::ZZ g0 { NTL::GCD(L_val, p02_1) };
        if (!NTL::divide(p0_zz-1, g0))
            continue;

        for (size_t j {i+1}; j < num_primes; ++j)
        {
            const long p1 { primes[j] };
            const NTL::ZZ p1_zz { primes[j] };
            if (NTL::divide(L_val, p1_zz)) 
                continue;

            const NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };
            if (NTL::divide(L_val, p12_1)) 
                continue;

            const NTL::ZZ g1 { NTL::GCD(L_val, p12_1) };
            if (!NTL::divide(p1_zz-1, g1))
                continue;

            if (NTL::GCD(p0_zz, p12_1) != 1 || NTL::GCD(p1_zz, p02_1) != 1)
                continue;

            const NTL::ZZ p0_inv { NTL::InvMod(p0_zz, p02_1) };
            const NTL::ZZ p1_inv { NTL::InvMod(p1_zz, p12_1) };
            const NTL::ZZ gcd { NTL::GCD(p02_1, p12_1) };

            NTL::ZZ r1, r2;
            NTL::rem(r1, p0_inv, gcd);
            NTL::rem(r2, p1_inv, gcd);
            if (r1 != r2)
                continue;

            factors.push_back(std::move(std::array<long, 2>{p0, p1}));
        }
    }
}



void
construct_primes_set(std::vector<long> &primes, const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &L_val, const Factorization &L, long max)
{
    NTL::ZZ p0_zz { nonrigid_factors[0] };
    NTL::ZZ p1_zz { nonrigid_factors[1] };
    NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };
    NTL::ZZ p_max {NTL::SqrRoot(L_val+1)};

#if LOG_LEVEL >= 1
    size_t count = 0;
    std::cout << "upper bound on prime <= " << p_max << "\n";
    if (max && max <= p_max)
        std::cout << "filtering primes <= " << max << " ...\n";
    else
        std::cout << "filtering primes <= " << p_max << " ...\n";
#endif

    NTL::PrimeSeq s;
    const size_t num_factors { L.primes.size() };
    size_t i { 0 };

    long p_prev { 0 };
    long p = s.next();
    while ( p != 0 && p <= p_max && (!max || p <= max))
    {
        /* check to make sure p is not a factor of L */
        for(; i < num_factors && L.primes[i] < p; ++i)
            ; /* NOP */

        if (i >= num_factors)
            break;
        else if (p != L.primes[i])
        {
            NTL::ZZ p_zz { p };
            if (NTL::GCD(p_zz, p02_1) == 1 && NTL::GCD(p_zz, p12_1) == 1)
            {
                NTL::ZZ p2_1 { NTL::sqr(p_zz)-1 };
                if(NTL::divide(L_val, p2_1))
                    primes.push_back(p);
            }
        }
        p_prev = p;
        p = s.next();
#if LOG_LEVEL >= 2
        if((count++ & STEP_MASK) == 0)
            std::cerr << "count: " << count << "\r";
#endif
    }

    /* NOTE: assming L does not contain a factor larger than what NTL::PrimeSeq supports */
    /* don't bother checking if p is a factor of L */
    while(p != 0 && p <= p_max && (!max || p <= max))
    {
        NTL::ZZ p_zz { p };
        if (NTL::GCD(p_zz, p02_1) == 1 && NTL::GCD(p_zz, p12_1) == 1)
        {
            if(NTL::divide(L_val, NTL::sqr(p_zz)-1))
            {
                primes.push_back(p);
            }
        }
        p_prev = p;
        p = s.next();
#if LOG_LEVEL >= 2
        if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
    }

#if LOG_LEVEL >= 2
    std::cerr << "count: " << count << "\n";
#endif

    if (p == 0 && p < p_max)
    {
        std::cout << "reached small prime limit, starting from p=" << p_prev << " ...\n";

        /* maxed out PrimeSeq, use NextPrime */
        NTL::ZZ p_zz {p_prev};
        NextPrime(p_zz, p_zz);
        while (p_zz < p_max && (!max || p_zz <= max))
        {
            /* check to make sure p is not a factor of L */
            for(; i < num_factors && L.primes[i] < p_zz; ++i); /* NOP */
            if (i >= num_factors) 
                break;
            else if(NTL::GCD(p_zz, p02_1) == 1 && NTL::GCD(p_zz, p12_1) == 1
                    && NTL::divide(L_val, NTL::sqr(p_zz)-1))
                primes.push_back(NTL::conv<long>(p_zz));

            NextPrime(p_zz, p_zz);
#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
        }

        /* don't bother checking if p_zz is a factor of L */
        while (p_zz < p_max && (!max || p_zz <= max))
        {
            if(NTL::GCD(p_zz, p02_1) == 1 && NTL::GCD(p_zz, p12_1) ==1
                    && NTL::divide(L_val, NTL::sqr(p_zz)-1))
                primes.push_back(NTL::conv<long>(p_zz));

            NextPrime(p_zz, p_zz);

#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
        }

#if LOG_LEVEL >= 2
        std::cerr << "count: " << count << "\n";
#endif
    }
}


// TODO: cache partial products
void
generate_a_values(std::vector<std::vector<long> > &a_values, const std::vector<long> &primes,
        const std::array<long,2> &nonrigid_factors, size_t min_terms, size_t max_terms)
{
    const long p0 { nonrigid_factors[0] };
    const long p1 { nonrigid_factors[1] };
    const NTL::ZZ p0_zz { p0 };
    const NTL::ZZ p1_zz { p1 };
    const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    const NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };

    NTL::ZZ_p::init(p12_1);
    NTL::ZZ_p p0_zzp {p0};
    std::cout << "1/p0 = " << NTL::inv(p0_zzp) << " (mod " << p12_1 << ")\n";

    NTL::ZZ_p::init(p02_1);
    NTL::ZZ_p p1_zzp {p1};
    std::cout << "1/p1 = " << NTL::inv(p1_zzp) << " (mod " << p02_1 << ")\n";

    size_t num_primes = primes.size();
    if (max_terms > 0)
        max_terms = MIN(num_primes, max_terms);
    else
        max_terms = num_primes;

    // go through all possible subset sizes starting from 2
    for (size_t t {MIN(min_terms,1)}; t <= max_terms; t++) 
    {
        size_t factors = t;
        std::vector<size_t> index_stack = {0};
#if LOG_LEVEL >= 1
        std::cout << "checking subsets of size " << t << " ...\n";
#endif

        // go through subsets of size factors suing a 'stack' (or odometer?)
        size_t index_size;
        while ((index_size = index_stack.size()) > 0) 
        {
            size_t top = index_size - 1;

            if (index_stack[top] < num_primes) 
            {
                if (index_size < factors) 
                {
                    index_stack.push_back(index_stack[top]+1);
                    continue;
                }

                // TODO: turn this into a loop
                NTL::ZZ_p::init(p02_1);
                NTL::ZZ_p prod { p1 };
                for(size_t i = 0; i < index_size; i++)
                    prod *= primes[index_stack[i]];

                if (NTL::IsOne(prod))
                {
                    NTL::ZZ_p::init(p12_1);
                    prod = p0;
                    for(size_t i = 0; i < index_size; i++)
                        prod *= primes[index_stack[i]];

                    //std::cout << "{";
                    //for(auto it=index_stack.cbegin(), end=index_stack.cend();
                    //        it != end; it++)
                    //    std::cout << primes[*it] << " ";
                    //std::cout << "} " << prod << "\n";
                    if(NTL::IsOne(prod))
                    {
                        std::vector<long> a_primes;
                        a_primes.reserve(t);
                        for(auto it=index_stack.cbegin(), end=index_stack.cend();
                                it != end; it++)
                            a_primes.push_back(primes[*it]);
                        std::cout << "found: " << a_primes << "\n";
                        a_values.push_back(std::move(a_primes));
                    }
                }

                index_stack[top]++;
            }
            else
            {
                index_stack.pop_back();
                if (top > 0) {
                    ++index_stack[top - 1];
                }

                continue;
            }
        }
    }
}

void
generate_cprimes(std::vector<std::vector<long>> &cprimes, const std::vector<long> primes, 
        const std::array<long, 2> &nonrigid_factors, 
        const NTL::ZZ &a_val, const std::vector<long> &a_primes,
        const NTL::ZZ &L_val, const Factorization &L,
        size_t min_terms, size_t max_terms)
{
    const long p0 { nonrigid_factors[0] };
    const long p1 { nonrigid_factors[1] };
    const NTL::ZZ p0_zz { p0 };
    const NTL::ZZ p1_zz { p1 };
    const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    const NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };
    const NTL::ZZ p0p1a { p0_zz * p1_zz * a_val };

    constexpr size_t num_bases { 3 };
    const std::array<NTL::ZZ, num_bases> prod_base {p02_1, p12_1, L_val};

    const size_t num_primes = primes.size();
    min_terms = MAX(min_terms, 2);
    if (max_terms > 0)
        max_terms = MIN(num_primes, max_terms);
    else
        max_terms = num_primes;

    // go through all possible subset sizes starting from 2
    for (size_t factors {min_terms}; factors <= max_terms; factors++) 
    {
        std::vector<size_t> index_stack = {0};
#if LOG_LEVEL >= 1
        std::cout << "checking subsets of size " << factors << " ...\n";
#endif

        /* cache of products */
        std::vector<std::array<NTL::ZZ, num_bases> > products; 
        products.reserve(factors);

        // go through subsets of size factors suing a 'stack' (or odometer?)
        size_t index_size;
        while ((index_size = index_stack.size()) > 0) 
        {
            const size_t prod_size  = products.size();
            size_t top = index_size - 1;
            size_t prod_top = prod_size - 1;

            if (index_stack[top] < num_primes) 
            {
                NTL::ZZ p_zz { primes[index_stack[top]] };
                if (prod_size+1 < factors) 
                {
                    std::array<NTL::ZZ, num_bases> prods;
                    for(size_t j { 0 }; j < num_bases; ++j)
                    {
                        NTL::ZZ prod { 1 };
                        if (prod_size)
                            prod = products[prod_top][j];
                        NTL::MulMod(prods[j], prod, p_zz, prod_base[j]);
                    }
                    products.push_back(std::move(prods));
                }

                if (index_size < factors) 
                {
                    index_stack.push_back(index_stack[top]+1);
                    continue;
                }
                /* IMPORTANT: must have at least TWO factors */
                /*            otherwise segfault  */
                std::array<NTL::ZZ, num_bases> &prods { products[top-1] };
                NTL::ZZ n0;

                // TODO: nested if's kinda ugly, refactor
                /* check n0 = 1 (mod p0^2-1) */
                NTL::MulMod(n0, prods[0], p_zz, prod_base[0]);
                if (NTL::IsOne(n0))
                {
                    /* check n0 = 1 (mod p1^2-1) */
                    NTL::MulMod(n0, prods[1], p_zz, prod_base[1]);
                    if (NTL::IsOne(n0))
                    {
                        /* check n0 * p0 * p1 * a = 1 (mod L) */
                        NTL::MulMod(n0, prods[1], p_zz, prod_base[2]);
                        NTL::MulMod(n0, n0, p0p1a, prod_base[2]);
                        if (NTL::IsOne(n0))
                        {
                            std::vector<long> n;
                            n.reserve(factors);
                            for(auto it=index_stack.cbegin(), end=index_stack.cend();
                                    it != end; ++it)
                                n.push_back(primes[*it]);
                            cprimes.push_back(std::move(n));
                            std::cout << "found: " << n << "\n";
                        }
                    }
                }

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
                        NTL::ZZ p_zz {primes[i]};
                        std::array<NTL::ZZ, num_bases> prods;
                        for(size_t j { 0 }; j < num_bases; ++j)
                        {
                            NTL::ZZ prod { 1 };
                            if (prod_size > 1)
                                prod = products[prod_top-1][j];
                            NTL::MulMod(prods[j], prod, p_zz, prod_base[j]);
                        }
                        products.push_back(std::move(prods));
                    }
                }

                continue;
            }
        }
    }
}
