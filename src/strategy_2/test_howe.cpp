#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <cmath>
#include <fstream>
#include <algorithm>

#include <NTL/ZZ.h>

#include <util.h>
#include <counting_factors.h>
#include <subset_product.h>
#include <strategy_2/nonrigid.h>

typedef std::unordered_map<std::array<NTL::ZZ, 2>, std::vector<std::vector<bool>>, ArrayHasher<NTL::ZZ, 2>> map_type_2;

void
test_construct_primes_set(std::vector<long> &primes, const long p0,
        const NTL::ZZ &L_val, const Factorization &L, long max)
{
    NTL::ZZ p0_zz { p0 };
    NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    NTL::ZZ p_max {NTL::SqrRoot(L_val+1)};
    if(max && max < p_max)
        p_max = max;

#if LOG_LEVEL >= 1
    size_t count = 0;
    std::cout << "upper bound on prime <= " << NTL::SqrRoot(L_val+1) << "\n";
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
            if (NTL::GCD(p_zz, p02_1) == 1)
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
        if (NTL::GCD(p_zz, p02_1) == 1)
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
        NextPrime(p_zz, p_zz+1);
        while (p_zz < p_max && (!max || p_zz <= max))
        {
            /* check to make sure p is not a factor of L */
            for(; i < num_factors && L.primes[i] < p_zz; ++i); /* NOP */
            if (i >= num_factors) 
                break;
            else if(NTL::GCD(p_zz, p02_1) == 1 && NTL::divide(L_val, NTL::sqr(p_zz)-1))
                primes.push_back(NTL::conv<long>(p_zz));

            NextPrime(p_zz, p_zz+1);
#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
        }

        /* don't bother checking if p_zz is a factor of L */
        while (p_zz < p_max && (!max || p_zz <= max))
        {
            if(NTL::GCD(p_zz, p02_1) == 1 && NTL::divide(L_val, NTL::sqr(p_zz)-1))
                primes.push_back(NTL::conv<long>(p_zz));

            NextPrime(p_zz, p_zz+1);

#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
        }

#if LOG_LEVEL >= 2
        std::cerr << "count: " << count << "\n";
#endif
    }
}


void
test_gen_cprimes_2way_all(
        const std::vector<long> &primes, 
        const long &nonrigid_factor,
        const NTL::ZZ &L_val,
        size_t min_size, size_t max_size)
{
    const long p0 { nonrigid_factor };
    const NTL::ZZ p0_zz { p0 };
    const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };

    constexpr size_t num_bases { 2 };
    std::array<NTL::ZZ, num_bases> prod_base {p02_1, L_val};
    map_type_2 map;

    /* split primes into two vectors */
    std::vector<long> h1_primes, h2_primes;
    split_half<long>(h1_primes, h2_primes, primes);
    const size_t h1_primes_size {h1_primes.size()};

    /* calculate group size */
    const NTL::ZZ lcm {get_lcm<std::array<NTL::ZZ, num_bases> >(prod_base, num_bases)};
    const NTL::ZZ mod_G {eulers_toitent(lcm)};

    /* clip min_size and max_size */
    bound<size_t>(min_size, 1, h1_primes_size);
    bound<size_t>(max_size, min_size, h1_primes_size);

    /* reserve space in hashmap */
    size_t num_subsets {calc_max_subsets(h1_primes_size, min_size, max_size)};
    constexpr size_t max_size_t {(size_t)-1};
    const size_t map_capacity {mod_G > max_size_t ? num_subsets 
        : MIN(num_subsets, NTL::conv<size_t>(mod_G)) };

#if LOG_LEVEL >= 1
    std::cout << "|G| = " << mod_G << " (" << ceil(NTL::log(mod_G)/log(2)) << " bits)\n";
    std::cout << "subset sizes " << "[" << min_size << "," << max_size << "]: "
        << num_subsets << " subsets\n";
    std::cout << "storing inverse subset products of first half (size " << h1_primes_size << ") " << " ...\n";
    std::cout << "map capacity: " << map_capacity << "\n";
    size_t inv_count {0};
#if LOG_LEVEL >= 2
    size_t count {0};
#endif
#endif
    map.reserve(map_capacity);
    /* go through every subset in first half, store inverse in hashmap */
    subsetprod_mod<2>(h1_primes, prod_base,
            [&](std::array<NTL::ZZ, num_bases>& prod_cache, 
                const std::vector<size_t> &indicies, size_t insert_index)->int
            {
                NTL::MulMod(prod_cache[1], prod_cache[1], p0_zz, prod_base[1]);
                std::array<NTL::ZZ, num_bases> invs;
                for(size_t i {0}; i < num_bases; ++i)
                   invs[i] = NTL::InvMod(prod_cache[i], prod_base[i]);

                // create vector<bool> representing this subset
                // and store in map
                std::vector<std::vector<bool>> &vec {map[invs]};

                std::vector<bool> subset;
                subset.reserve(h1_primes_size);
                subset.resize(h1_primes_size, false);
                for (auto i : indicies) subset[i] = true;

                vec.push_back(std::move(subset));
#if LOG_LEVEL >= 1
                if(vec.size() == 1)
                    inv_count++;
#if LOG_LEVEL >= 2
                if((count++ & STEP_MASK) == 0)
                    std::cerr << "count: " << count << "\r";
#endif
#endif
                return 1;
            }, min_size, max_size);

#if LOG_LEVEL >= 1
#if LOG_LEVEL >= 2
    std::cerr << "count: " << count << "\r";
    count = 0;
#endif
    std::cout << "found " << inv_count << " distinct subset products\n";
    std::cout << "checking for inverse subset products of second half (size " << h2_primes.size() << ") ...\n";
#endif
    /* go through second half, check whether inverses exist in the hashmap */
    subsetprod_mod<2>(h2_primes, prod_base,
            [&](std::array<NTL::ZZ, num_bases>& prod_cache, 
                const std::vector<size_t> &indicies, size_t insert_index)->int
            {
                auto it {map.find(prod_cache)};
                if (it != map.end())
                {
                    for (const std::vector<bool> &other_half : it->second)
                    {
                        std::vector<long> cprimes;
                        const size_t other_half_size {other_half.size()};

                        // assert other_half.size() == h1_primes.size()
                        for (size_t i = 0; i < other_half_size; i++)
                            if (other_half[i])
                                cprimes.push_back(h1_primes[i]);

                        for(auto i : indicies)   cprimes.push_back(h2_primes[i]);
                        std::cout << "found: " << cprimes << "\n";
                    }
                }
#if LOG_LEVEL >= 2
                if((count++ & STEP_MASK) == 0)
                    std::cerr << "count: " << count << "\r";
#endif
                return 0;
            }, min_size, max_size);

#if LOG_LEVEL >= 2
    std::cerr << "count: " << count << "\r";
#endif
}


/**
 * Format for commandline arguments(<> required, [] optional):
 *  [1]     <max> <p0> <min_terms> <max_terms>
 */
int
main(int argc, char **argv)
{
    const long max { NTL::conv<long>(argv[1]) };
    const long p0 { 1153 };
    const NTL::ZZ p0_zz {p0};
    const NTL::ZZ p02_1 {NTL::sqr(p0_zz)-1};

    const size_t min_terms { NTL::conv<size_t>(argv[2]) };
    const size_t max_terms { NTL::conv<size_t>(argv[3]) };

    Factorization L { parse_factorization("2^7 3^3 5^2 7^1 11^1 13^1 17^1 19^1 29^1 31^1") };

    NTL::ZZ L_val;
    multiply_factors(L_val, L.primes, L.powers);
    std::cout << "L = " << L_val << "\n";

    std::cout << "p0 = " << p0 << "\n";
    std::cout << "p0^2-1 = " << p02_1 << "\n";

    std::vector<long> primes_set;
    test_construct_primes_set(primes_set, p0, L_val, L, max);
    std::cout << "P = " << primes_set << "\n";

    const size_t primes_set_size { primes_set.size() };
    std::cout << "|P| = " << primes_set_size << "\n";

    test_gen_cprimes_2way_all(primes_set, p0, L_val, min_terms, max_terms);
}
