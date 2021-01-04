#include <iostream>
#include <vector>
#include <functional>
#include <unordered_map>
#include <map>
#include <climits>
#include <cstdlib>
#include <cmath>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <config.h>
#include <subset_product.h>
#include <counting_factors.h>
#include <util.h>
#include <strategy_2/nonrigid.h>


template <size_t N>
using map_type = std::unordered_map<std::array<NTL::ZZ, N>, std::vector<std::vector<bool>>, ArrayHasher<NTL::ZZ, N>>;


/**
 * Generate possible choices for non-rigid factors for the given choice of L.
 * The pair of non-rigid factors stored as array<long,2> and added to
 * `factors` vector.
 *
 * @param   factors             vector to add the generate non-rigid factor pairs
 * @param   L_val               value of parameter L as NTL::ZZ
 * @param   max                 maximum integer to consider
 */
void
generate_possible_factors(std::vector<std::array<long,2> > &factors, const NTL::ZZ &L_val, const long max)
{
    std::vector<long> primes;
#if LOG_LEVEL >= 1
    std::cout << "generating primes up to " << max << " ...\n";
#if LOG_LEVEL >= 2
    size_t count {0};
    size_t cond[4] {0, 0, 0, 0};

#endif
#endif

    NTL::ZZ p { 1 };
    while (p < max)
    {
        NTL::NextPrime(p, p+1);
#if LOG_LEVEL >= 2
        if ((count++ & STEP_MASK) == 0)
            std::cerr << " count: " << count
                << "; cond: " << cond[0] << ", " << cond[1] << "\r";
#endif

        if (NTL::divide(L_val, p))
            continue;

        const NTL::ZZ p2_1 { NTL::sqr(p)-1 };
        if (NTL::divide(L_val, p2_1))
            continue;
#if LOG_LEVEL >= 2
        cond[0]++;
#endif

        const NTL::ZZ g { NTL::GCD(L_val, p2_1) };
        if (!NTL::divide(p-1, g))
            continue;
#if LOG_LEVEL >= 2
        cond[1]++;
#endif

        primes.push_back(NTL::conv<long>(p));
    }

#if LOG_LEVEL >= 1
#if LOG_LEVEL >= 2
    std::cerr << " count: " << count
        << "; cond: " << cond[0] << ", " << cond[1] << "\n";
    count = 0;
#endif
    std::cout << "filtering pairs of primes from list of size " << primes.size() << " ...\n";
#endif

    const size_t num_primes { primes.size() };
    for(size_t i { 0 }; i < num_primes; ++i)
    {
        const long p0 { primes[i] };
        const NTL::ZZ p0_zz { primes[i] };
        const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };

        for (size_t j {i+1}; j < num_primes; ++j)
        {
            const long p1 { primes[j] };
            const NTL::ZZ p1_zz { primes[j] };
            const NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };

#if LOG_LEVEL >= 2
            if ((count++ & STEP_MASK) == 0)
                std::cerr << " count: " << count
                    << "; cond: " << cond[2] << ", " << cond[3] << "\r";
#endif

            if (NTL::GCD(p0_zz, p12_1) != 1 || NTL::GCD(p1_zz, p02_1) != 1)
                continue;
#if LOG_LEVEL >= 2
            cond[2]++;
#endif

            const NTL::ZZ p0_inv { NTL::InvMod(p0_zz, p02_1) };
            const NTL::ZZ p1_inv { NTL::InvMod(p1_zz, p12_1) };
            const NTL::ZZ gcd { NTL::GCD(p02_1, p12_1) };

            NTL::ZZ r1, r2;
            NTL::rem(r1, p0_inv, gcd);
            NTL::rem(r2, p1_inv, gcd);
            if (r1 != r2)
                continue;
#if LOG_LEVEL >= 2
            cond[3]++;
#endif

            std::cout << "found: " << p0 << " " << p1 << "\n";
            factors.push_back(std::move(std::array<long, 2>{p0, p1}));
        }
    }
#if LOG_LEVEL >= 2
    std::cerr << " count: " << count
        << "; cond: " << cond[2] << ", " << cond[3] << "\n";
#endif
}


/**
 * Generate the primes in P(2,L) given the choice of L, p0, and p1.
 *
 * @param   primes              vector to add the generate primes
 * @param   nonrigid_factors    array<NTL::ZZ,2> of non-rigid factors {p0,p1}
 * @param   L_val               value of parameter L as NTL::ZZ
 * @param   L                   the factorization of L as a `Factorization` object
 * @param   max                 maximum integer to consider
 */
void
construct_primes_set(std::vector<long> &primes, const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &L_val, const Factorization &L, long max)
{
    NTL::ZZ p0_zz { nonrigid_factors[0] };
    NTL::ZZ p1_zz { nonrigid_factors[1] };
    NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };
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
    while ( p != 0 && p <= p_max)
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
    while(p != 0 && p <= p_max)
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
        NextPrime(p_zz, p_zz+1);
        while (p_zz < p_max)
        {
            /* check to make sure p is not a factor of L */
            for(; i < num_factors && L.primes[i] < p_zz; ++i); /* NOP */
            if (i >= num_factors)
                break;
            else if(NTL::GCD(p_zz, p02_1) == 1 && NTL::GCD(p_zz, p12_1) == 1
                    && NTL::divide(L_val, NTL::sqr(p_zz)-1))
                primes.push_back(NTL::conv<long>(p_zz));

            NextPrime(p_zz, p_zz+1);
#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
        }

        /* don't bother checking if p_zz is a factor of L */
        while (p_zz < p_max)
        {
            if(NTL::GCD(p_zz, p02_1) == 1 && NTL::GCD(p_zz, p12_1) ==1
                    && NTL::divide(L_val, NTL::sqr(p_zz)-1))
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


/**
 * Go through subsets of `primes` to look for possible choices for paramter a
 * and store them in `a_values` given the choice of L, p0, and p1.
 * Uses 2-set subset-product algorithm.
 *
 * @param   a_values            vector to store possible choices for parameter a
 * @param   primes              set of primes, P(2,L)
 * @param   nonrigid_factors    array<NTL::ZZ,2> of non-rigid factors {p0,p1}
 * @param   L_val               value of parameter L as NTL::ZZ
 * @param   min_size            minium size of subset to consider for each of
 *                              the two paritions of P(2,L)
 * @param   min_size            maximum size of subset to consider for each of
 *                              the two paritions of P(2,L)
 */
void
generate_a_values(std::vector<std::vector<long> > &a_values, const std::vector<long> &primes,
        const std::array<long,2> &nonrigid_factors,
        const NTL::ZZ &L_val,
        size_t min_terms, size_t max_terms)
{
    const long p0 { nonrigid_factors[0] };
    const long p1 { nonrigid_factors[1] };
    const NTL::ZZ p0_zz { p0 }, p1_zz {p1};
    const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    const NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };
    const NTL::ZZ gcd { NTL::GCD(NTL::GCD(p02_1, p12_1), L_val)};
    const NTL::ZZ p0p1 {NTL::MulMod(p0_zz, p1_zz, gcd)};

    constexpr size_t num_bases {2};
    const std::array<NTL::ZZ, num_bases> prod_base {p02_1, p12_1};
    map_type<2> map;

    /* split primes into two vectors */
    std::vector<long> h1_primes, h2_primes;
    split_half(h1_primes, h2_primes, primes);
    const size_t h1_primes_size {h1_primes.size()};
    const size_t num_primes {primes.size()};

#if LOG_LEVEL >= 1
    const NTL::ZZ inv_p0 {NTL::InvMod(p0_zz, p12_1)};
    const NTL::ZZ inv_p1 {NTL::InvMod(p1_zz, p02_1)};
    std::cout << "\np0^-1 = " << inv_p0 << " (mod " << p12_1 << ")\n";
    std::cout << "p1^-1 = " << inv_p1 << " (mod " << p02_1 << ")\n";
    std::cout << "gcd(p0^2-1,p1^2-1,L) = " << gcd << "\n";

    std::cout << "storing inverse subset products of first half (size " << h1_primes_size << ") " << " ...\n";
    size_t inv_count {0};
#if LOG_LEVEL >= 2
    size_t count {0};
#endif
#endif

    /* go through every subset of first half, store inverses * target in hashmap */
    subsetprod_mod<2>(h1_primes, prod_base,
            [&](std::array<NTL::ZZ, num_bases>& prod_cache,
                const std::vector<size_t> &indicies, size_t insert_index)->int
            {
                NTL::MulMod(prod_cache[0], prod_cache[0], p1_zz, prod_base[0]);
                NTL::MulMod(prod_cache[1], prod_cache[1], p0_zz, prod_base[1]);

                /* calculate inverses */
                std::array<NTL::ZZ, num_bases> invs;
                for(size_t j{0}; j < num_bases; ++j)
                    NTL::InvMod(invs[j], prod_cache[j], prod_base[j]);

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
                    std::cerr << "count: " << count
                        << ", inv: " << inv_count << "\r";
#endif
#endif
                return 1;
            }, min_terms, max_terms);

#if LOG_LEVEL >= 1
#if LOG_LEVEL >= 2
        std::cerr << "count: " << count << "\r";
        count = 0;
#endif
    std::cout << "found " << inv_count << " distinct subset products\n";
    std::cout << "checking for inverse subset products of second half (size " << h2_primes.size() << ") ...\n";
#endif

    subsetprod_mod<2>(h2_primes, prod_base,
            [&](std::array<NTL::ZZ, num_bases>& prod_cache,
                const std::vector<size_t> &indicies, size_t insert_index)->int
            {
                const std::vector<std::vector<bool>> &ret {map[std::move(prod_cache)]};
                if (ret.size() > 0)
                {
                    for (const std::vector<bool> &other_half: ret)
                    {
                        std::vector<long> a_primes;
                        const size_t other_half_size {other_half.size()};
                        a_primes.reserve(other_half_size + indicies.size());

                        // assert other_half.size() == h1_primes.size()
                        for (size_t i = 0; i < other_half_size; i++)
                            if (other_half[i])
                                a_primes.push_back(h1_primes[i]);
                        for(auto i : indicies)   a_primes.push_back(h2_primes[i]);

                        /* check that p0p1a^-1 = 1 (mod gcd(p0^2-1, p1^2-1, L)) */
                        NTL::ZZ p0p1a {p0p1};
                        for (auto v : a_primes)
                            NTL::MulMod(p0p1a, p0p1a, v, gcd);
                        if(!NTL::IsOne(NTL::InvMod(p0p1a, gcd)))
                            continue;

                        std::cout << a_primes.size() << ",|P| = " << num_primes-a_primes.size()
                            << "," << a_primes << "\n";
                        a_values.push_back(std::move(a_primes));
                    }
                }
#if LOG_LEVEL >= 2
                if((count++ & STEP_MASK) == 0)
                    std::cerr << "count: " << count << "\r";
#endif
                return 0;
            }, min_terms, max_terms);
}


NTL::ZZ&
calc_target_residue(NTL::ZZ &target, const NTL::ZZ &L_val,
        const NTL::ZZ &a_val,
        const std::array<long,2> &nonrigid_factors)
{
    const NTL::ZZ p0_zz { nonrigid_factors[0] };
    const NTL::ZZ p1_zz { nonrigid_factors[1] };
    const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    const NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };
    const NTL::ZZ p0p1a { p0_zz * p1_zz * a_val };

    const NTL::ZZ lcm_p0p1 {
        get_lcm<std::array<NTL::ZZ,2>>(std::array<NTL::ZZ,2>{p02_1, p12_1},2)
    };
    // std::cout << "lcm(p0^2-1,  p1^2-1) = " << lcm_p0p1 << "\n";

    {
        const NTL::ZZ gcd { NTL::GCD(L_val, lcm_p0p1) };
        const NTL::ZZ base { lcm_p0p1/gcd };

        NTL::ZZ Ld_inv;
        NTL::rem(Ld_inv, L_val/gcd, base);
        NTL::InvMod(Ld_inv, Ld_inv, base);
        // std::cout << "(L/gcd)^-1 = " << Ld_inv << " (mod " << base << ")\n";

        NTL::ZZ p0p1a_inv;
        NTL::rem(p0p1a_inv, p0p1a, base);
        NTL::InvMod(p0p1a_inv, p0p1a_inv, base);
        // std::cout << "(p0p1a)^-1 = " << p0p1a << " (mod " << base << ")\n";

        NTL::ZZ t {1-p0p1a_inv};
        if (t < 0) t = base - t;
        t = (t / gcd) * Ld_inv;
        target = t;
    }

    std::cout << "t1 = " << target << "\n";
    return target;
}


/**
 * Generate order2-carmichael numbers using 2-set subset-product algorithm.
 * of `primes`.
 * @param   primes              set of primes, P(2,L)
 * @param   nonrigid_factors    array<NTL::ZZ,2> of non-rigid factors {p0,p1}
 * @param   a_val               parameter a for the given choice of L, p0, and p1
 * @param   L_val               value of parameter L as NTL::ZZ
 * @param   min_size            minium size of subset to consider for each of
 *                              the two paritions of P(2,L)
 * @param   min_size            maximum size of subset to consider for each of
 *                              the two paritions of P(2,L)
 */
void
gen_cprimes_2way_all(
        const std::vector<long> &primes,
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        size_t min_size, size_t max_size)
{
    const long p0 { nonrigid_factors[0] };
    const long p1 { nonrigid_factors[1] };
    const NTL::ZZ p0_zz { p0 };
    const NTL::ZZ p1_zz { p1 };
    const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    const NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };
    const NTL::ZZ p0p1a { p0_zz * p1_zz * a_val };

    constexpr size_t num_bases { 3 };
    std::array<NTL::ZZ, num_bases> prod_base {p02_1, p12_1, L_val};
    map_type<3> map;

    /* split primes into two vectors */
    std::vector<long> h1_primes, h2_primes;
    split_half<long>(h1_primes, h2_primes, primes);
    const size_t h1_primes_size {h1_primes.size()};

    /* calculate group size */
    const NTL::ZZ lcm {get_lcm<std::array<NTL::ZZ, num_bases> >(prod_base, num_bases)};
    const NTL::ZZ mod_G {eulers_toitent(lcm)};

    /* clip min_size and max_size */
    min_size = bound<size_t>(min_size, 1, h1_primes_size);
    max_size = bound<size_t>(max_size, min_size, h1_primes_size);

    /* reserve space in hashmap */
    size_t num_subsets {calc_max_subsets(h1_primes_size, min_size, max_size)};
    constexpr size_t max_size_t {(size_t)-1};
    const size_t map_capacity {mod_G > max_size_t ? num_subsets
        : MIN(num_subsets, NTL::conv<size_t>(mod_G)) };

#if LOG_LEVEL >= 1
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
    subsetprod_mod<3>(h1_primes, prod_base,
            [&](std::array<NTL::ZZ, num_bases>& prod_cache,
                const std::vector<size_t> &indicies, size_t insert_index)->int
            {
                NTL::MulMod(prod_cache[2], prod_cache[2], p0p1a, prod_base[2]);
                std::array<NTL::ZZ, num_bases> invs;
                for(size_t i {0}; i < num_bases; ++i)
                   invs[i] = NTL::InvMod(prod_cache[i], prod_base[i]);

                // create vector<bool> representing this subset
                // and store in map
                std::vector<std::vector<bool>> &vec {map[invs]};

                std::vector<bool> subset;
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
    subsetprod_mod<3>(h2_primes, prod_base,
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
 * Generate order2-carmichael numbers by picking random subsets
 * of `primes`.
 * @param   primes              set of primes, P(2,L)
 * @param   nonrigid_factors    array<NTL::ZZ,2> of non-rigid factors {p0,p1}
 * @param   a_val               parameter a for the given choice of L, p0, and p1
 * @param   L_val               value of parameter L as NTL::ZZ
 * @param   num_trials          number times to pick random subsets.
 *                              If zero or > |G|, then |G| is chosen.
 * @param   min_size            minium size of subset to generate
 * @param   max_size            maximum size of subset to generate
 */
void
gen_cprimes_2way_prob(
        const std::vector<long> &primes,
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        size_t num_trials,
        size_t min_size, size_t max_size)
{
    const long p0 { nonrigid_factors[0] };
    const long p1 { nonrigid_factors[1] };
    const NTL::ZZ one {1};
    const NTL::ZZ p0_zz { p0 };
    const NTL::ZZ p1_zz { p1 };
    const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    const NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };
    const NTL::ZZ p0p1a { p0_zz * p1_zz * a_val };

    constexpr size_t num_bases { 3 };
    std::array<NTL::ZZ, num_bases> prod_base {p02_1, p12_1, L_val};
    std::unordered_map<std::array<NTL::ZZ, 3>, std::vector<size_t>, ArrayHasher<NTL::ZZ, 3>> map;

    /* calculate the total group size */
    const NTL::ZZ mag_G {eulers_toitent(get_lcm<std::array<NTL::ZZ, num_bases> >(prod_base, num_bases))};

    /* split primes into two vectors */
    std::vector<long> h1_primes, h2_primes;
    split_half(h1_primes, h2_primes, primes);
    const size_t h1_primes_size {h1_primes.size()};
    const size_t h2_primes_size {h2_primes.size()};

    const size_t sqrt_g {NTL::conv<size_t>(NTL::SqrRoot(mag_G))};
    const size_t max_num_subsets {calc_max_subsets(h1_primes_size, min_size, max_size)};
    if (!num_trials)
    {
        num_trials  = MIN(sqrt_g, max_num_subsets);
    }

#if LOG_LEVEL >= 2
    size_t count {0};
#if LOG_LEVEL >= 1
    std::cout << "|G| = " << mag_G << " (" << ceil(NTL::log(mag_G)/log(2)) << " bits)\n";
    std::cout << "randomly choosing " << num_trials << " subsets from first half (size "
        << h1_primes_size << ")...\n";
    size_t inv_count {0};
#endif
#endif

    size_t subset_count {0};
    std::vector<std::vector<size_t> > subsets;
    subsets.reserve(num_trials);

    for(size_t i {0}; i < num_trials; ++i)
    {
        std::unique_ptr<std::vector<size_t>> indicies_ptr {random_subset(h1_primes_size)};
        std::vector<size_t> &indicies { *indicies_ptr };

        /* calculate the product of the random subset */
        NTL::ZZ prods[num_bases];
        for(size_t j {0}; j < num_bases; ++j)
            NTL::MulMod(prods[j], one, h1_primes[indicies[0]], prod_base[j]);
        for(size_t m {1}, len=indicies.size(); m < len; ++m)
        {
            for(size_t j {0}; j < num_bases; ++j)
                NTL::MulMod(prods[j], prods[j], h1_primes[indicies[m]], prod_base[j]);
        }
        NTL::MulMod(prods[2], prods[2], p0p1a, prod_base[2]);


        /* find the inverses */
        std::array<NTL::ZZ, num_bases> invs;
        for(size_t j {0}; j < num_bases; ++j)
            NTL::InvMod(invs[j], prods[j], prod_base[j]);

        /* add inverse to add with index as value */
        std::vector<size_t> &vec {map[std::move(invs)]};
        vec.push_back(subset_count);

        /* add random subset to subsets list */
        subsets.push_back(std::move(indicies));
        subset_count++;
#if LOG_LEVEL >= 2
        if((count++ & STEP_MASK) == 0)
            std::cerr << "count: " << count << "\r";
#if LOG_LEVEL >= 1
        if (vec.size() == 1)
            inv_count++;
#endif
#endif
    }

#if LOG_LEVEL >= 2
    std::cerr << "count: " << subset_count << "\r";
    count = 0;
#if LOG_LEVEL >= 1
    std::cout << "found " << inv_count << " distinct subset products\n";
    std::cout << "randomly choosing " << num_trials << " subsets from second half ...\n";
#endif
#endif

    for(size_t i {0}; i < num_trials; ++i)
    {
        std::unique_ptr<std::vector<size_t>> indicies_ptr {random_subset(h2_primes_size)};
        std::vector<size_t> &indicies { *indicies_ptr };

        /* calculate the product of the random subset */
        std::array<NTL::ZZ, num_bases> prods;
        for(size_t j {0}; j < num_bases; ++j)
            NTL::MulMod(prods[j], one, h1_primes[indicies[0]], prod_base[j]);
        for(size_t m {1}, len=indicies.size(); m < len; ++m)
        {
            for(size_t j {0}; j < num_bases; ++j)
                NTL::MulMod(prods[j], prods[j], h1_primes[indicies[m]], prod_base[j]);
        }

        /* check if match exists in map */
        const std::vector<size_t> &vec {map[std::move(prods)]};
        if (vec.size() > 0)
        {
            for (auto ind : vec)
            {
                std::vector<long> cprimes;
                const std::vector<size_t> &other_half {subsets[ind]};
                for(auto i : other_half) cprimes.push_back(h1_primes[i]);
                for(auto i : indicies)   cprimes.push_back(h2_primes[i]);
                std::cout << "found: " << cprimes << "\n";
            }
        }
#if LOG_LEVEL >= 2
        if((count++ & STEP_MASK) == 0)
            std::cerr << "count: " << count << "\r";
#endif
    }

#if LOG_LEVEL >= 2
    std::cerr << "count: " << count << "\r";
#endif
}


/*
 * Generate order2-carmichael numbers with two non-rigid factors given in
 * `nonrigid_factors` using 4-set algorithm.
 *
 * @param  primes               set of primes P(2,L)
 * @param  nonrigid_factors     array<NTL:ZZ,2> of non-rigid factors {p0, p1}
 * @param  a_val                parameter a for given choice for L, p0, and p1
 * @param  L_val                value for paramter L as NTL::ZZ
 * @param  min_size             minimum size of subsets to consider for each of
 *                              of the four partitions of `primes`
 * @param  min_size             maximum size of subsets to consider for each of
 *                              of the four partitions of `primes`
 */
void
gen_cprimes_4way_all(
        const std::vector<long> &primes,
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        size_t min_size, size_t max_size
    )
{
    const long p0 { nonrigid_factors[0] };
    const long p1 { nonrigid_factors[1] };
    const NTL::ZZ p0_zz { p0 };
    const NTL::ZZ p1_zz { p1 };
    const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    const NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };
    const NTL::ZZ p0p1a { p0_zz * p1_zz * a_val };

    /* split primes set into 4 */
    std::vector<std::vector<long>> primes_n {split_vector(primes, 4)};

    /* store the subset-product mod p02_1 and p12_1 as array<NTL::ZZ,2>
     * of every subset of each of the primes set in primes_n */
    std::array<NTL::ZZ, 2> prod_base {p02_1, p12_1};

    /* subset-products will be stored in a map with an array of
     * subset-products as the key and the subsets as vector<bool>
     * as the value */
    std::array<map_type<2>, 4> maps;
    for (size_t i = 0; i < 4; i++)
    {
        const auto &set = primes_n[i];
        const size_t set_size = set.size();
        const size_t pmin_size = bound<size_t>(min_size, 1, set_size);
        const size_t pmax_size = bound<size_t>(max_size, pmin_size, set_size);
        size_t num_subsets {calc_max_subsets(set.size(), pmin_size, pmax_size)};

#if LOG_LEVEL >= 1
        printf("(%lu) subset sizes [%lu,%lu] of %lu: %lu subsets\n", i, pmin_size,
                pmax_size, set_size, num_subsets);
        size_t inv_count {0};
#if LOG_LEVEL >= 2
        size_t count {0};
#endif
#endif
        maps[i].reserve(num_subsets);
        subsetprod_mod<2>(set, prod_base,
                [&](std::array<NTL::ZZ, 2> &prod_cache,
                    const std::vector<size_t> &indicies, size_t insert_index)->int
                {
                    /* store the inverse subset-product instead of subset-product
                     * mod prod_base for one of the primes sets in the pair */
                    if (i % 2  == 0)
                    {
                        /* primes_n[0] and primes_n[1] */
                        for (size_t j {0}; j < 2; ++j)
                            prod_cache[j] = NTL::InvMod(prod_cache[j], prod_base[j]);
                    }

                    /* create a vector<bool> representing this subset and add it to map */
                    std::vector<std::vector<bool>> &vec {maps[i][prod_cache]};

                    /* `indicies` indexes elements of `set` or `primes_n[i]` */
                    std::vector<bool> subset;
                    subset.resize(set_size, false);
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
                },
            pmin_size, pmax_size
        );

#if LOG_LEVEL >= 1
#if LOG_LEVEL >= 2
        printf("count: %lu\n", count);
#endif
        printf("inverses: %lu\n", inv_count);
#endif
    }

#if LOG_LEVEL >= 1
    printf("\njoining [0] and [1] ...\n");
    size_t count {0};
#endif
    /* next we care about subset-product mod L_val only */
    std::array<NTL::ZZ, 1> prod_base_2 {L_val};
    /* join subsets in maps[0] and map[1] that are in unitary subgroup mod {p02_1, p12_1}
     * and store them in map2_0 using their product mod L_val as the keys */
    map_type<1> map2_0;
    join<2,1>(map2_0, maps[0], maps[1], primes_n[0], primes_n[1], prod_base_2,
            [&](const std::array<NTL::ZZ,1> &prod0, const std::array<NTL::ZZ,1> &prod1)
                ->std::array<NTL::ZZ,1>
            {
                /* prod0 and prod1 are the product of each subset from map[0] and map[1]
                 * mod L_val; we return what they key for joined subset should be. */
                std::array<NTL::ZZ, 1> prod {NTL::ZZ{1}};
                /* In this case, they key is the product of the joined subset mod L_val */
                NTL::MulMod(prod[0], prod0[0], prod1[0], L_val);
#if LOG_LEVEL >= 1
                count++;
#endif
                return prod;
            }
        );
    /* the subsets we care about have been copied already, so clear them */
    maps[0].clear();
    maps[1].clear();

#if LOG_LEVEL >= 1
    printf("subset count: %lu\n", count);
    printf("\njoining [2] and [3] ...\n");
    count = 0;
#endif

    /* join subsets in map[2] and map[3] this time*/
    map_type<1> map2_1;
    join<2,1>(map2_1, maps[2], maps[3], primes_n[2], primes_n[3],prod_base_2,
            [&](const std::array<NTL::ZZ,1> &prod0, const std::array<NTL::ZZ,1> &prod1)
                ->std::array<NTL::ZZ,1>
            {
                /* this time they key for the joined subset is
                 * the inverset of (the product * p0 * p1 * a) mod L_val */
                std::array<NTL::ZZ, 1> prod {p0p1a};
                NTL::MulMod(prod[0], prod0[0], prod1[0], L_val);
                NTL::InvMod(prod[0], prod[0], L_val);
#if LOG_LEVEL >= 1
                count++;
#endif
                return prod;
            }
        );
    maps[2].clear();
    maps[3].clear();
#if LOG_LEVEL >= 1
    printf("subset count: %lu\n", count);
    printf("\nfinal join ...\n");
    count = 0;
#endif

    /* the subsets stored in maps2_0 and maps2_1 are stored as vector<bool> and
     * draw halfs of the original `primes` set */
    std::vector<long> primes01 {join_vector(primes_n[0], primes_n[1])};
    std::vector<long> primes23 {join_vector(primes_n[2], primes_n[3])};

    std::array<NTL::ZZ,0> empty_base;
    map_type<0> map_final;
    join<1,0>(map_final, map2_0, map2_1, primes01, primes23, empty_base,
            [&](const std::array<NTL::ZZ, 0> &prod0, const std::array<NTL::ZZ, 0> &prod1)
                ->std::array<NTL::ZZ,0>
            {
#if LOG_LEVEL >= 1
                count++;
#endif
                /* we don't care about the key in the final join, just store them under the same key */
                return empty_base;
            }
        );

#if LOG_LEVEL >= 1
    printf("subset count: %lu\n", count);
#endif
}



/*
 * Generate order2-carmichael numbers with two non-rigid factors given in
 * `nonrigid_factors` using 8-set algorithm.
 *
 * @param  primes               set of primes P(2,L)
 * @param  nonrigid_factors     array<NTL:ZZ,2> of non-rigid factors {p0, p1}
 * @param  a_val                parameter a for given choice for L, p0, and p1
 * @param  L_val                value for paramter L as NTL::ZZ
 * @param  min_size             minimum size of subsets to consider for each of
 *                              of the 8 partitions of `primes`
 * @param  min_size             maximum size of subsets to consider for each of
 *                              of the 8 partitions of `primes`
 */
void
gen_cprimes_8way_all(
        const std::vector<long> &primes,
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        size_t min_size, size_t max_size
    )
{
    const long p0 { nonrigid_factors[0] };
    const long p1 { nonrigid_factors[1] };
    const NTL::ZZ p0_zz { p0 };
    const NTL::ZZ p1_zz { p1 };
    const NTL::ZZ p02_1 { NTL::sqr(p0_zz)-1 };
    const NTL::ZZ p12_1 { NTL::sqr(p1_zz)-1 };
    const NTL::ZZ p0p1a { p0_zz * p1_zz * a_val };

    /* split primes set into 8 */
    std::vector<std::vector<long>> partitions {split_vector(primes, 8)};
    std::array<NTL::ZZ, 4> prod_bases {p02_1, p12_1, L_val, NTL::ZZ{1}};

    std::vector<map_type<1>> maps0;
    maps0.resize(8);
    for (size_t n=0; n < 8; n++)
    {
        const std::vector<long> &Pn { partitions[n] };
        const size_t Pn_size { Pn.size() };
        const size_t pmin_size = bound<size_t>(min_size, 1, Pn_size);
        const size_t pmax_size = bound<size_t>(max_size, pmin_size, Pn_size);
        size_t num_subsets {calc_max_subsets(Pn_size, pmin_size, pmax_size)};

#if LOG_LEVEL >= 1
        printf("(%lu) subset sizes [%lu,%lu] of %lu: %lu subsets\n", n, pmin_size,
                pmax_size, Pn_size, num_subsets);
        size_t inv_count {0};
#if LOG_LEVEL >= 2
        size_t count {0};
#endif
#endif

        maps0[n].reserve(num_subsets);
        /* calculate & store subset-product mod p0^2-1 for subsets of
         * sizes pmin_size to pmax_size for all of the 8 partitions */
        subsetprod_mod<1>(Pn, std::array<NTL::ZZ,1>{prod_bases[0]},
                [&](std::array<NTL::ZZ, 1> &prod_cache,
                    const std::vector<size_t> &indicies, size_t insert_index)->int
                {
                    /* store the inverse subset-product instead of subset-product
                     * mod p02_1 for one of the primes sets in the pair */
                    if (n % 2  == 0)
                    {
                        prod_cache[0] = NTL::InvMod(prod_cache[0], prod_bases[0]);
                    }

                    /* create a vector<bool> representing this subset and add it to map */
                    std::vector<std::vector<bool>> &vec {maps0[n][prod_cache]};

                    std::vector<bool> subset;
                    subset.resize(Pn_size, false);
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
                },
            pmin_size, pmax_size
        );
#if LOG_LEVEL >= 1
#if LOG_LEVEL >= 2
        fprintf(stderr, "count: %lu\r", count);
#endif

        fprintf(stderr, "inverses: %lu\n", inv_count);
#endif
    }

    std::vector<map_type<1>> maps1;
    std::vector<std::vector<long>> joined_partitions = partitions;
    for (size_t layer=1; layer <= 3; layer++)
    {
        std::vector<std::vector<long>> next_partitions;
        std::vector<map_type<1>> &prev_maps { layer % 2 == 1 ? maps0 : maps1};
        std::vector<map_type<1>> &next_maps { layer % 2 == 1 ? maps1 : maps0};
        next_maps.clear();

        const size_t num_joined {(size_t)1 << layer};
        const size_t num_joins  {(size_t)8 >> layer};
        next_maps.resize(num_joins);
        next_partitions.reserve(num_joins);

        for (size_t k=0; k < num_joins; k++)
        {
            const size_t k0 = 2*k; 
            const size_t k1 = k0+1;
#if LOG_LEVEL >= 1
            printf("\nlayer %lu, joining [%lu] and [%lu] ...\n", layer-1, k0, k1);
#endif
            join<1,1>(
                    next_maps[k],
                    prev_maps[k0], prev_maps[k1], 
                    joined_partitions[k0], joined_partitions[k1],
                    std::array<NTL::ZZ,1>{prod_bases[layer]},
                    [&](const std::array<NTL::ZZ,1> &prod0, const std::array<NTL::ZZ,1> &prod1)
                        ->std::array<NTL::ZZ,1>
                    {
                        std::array<NTL::ZZ, 1> prod {NTL::ZZ{1}};
                        if (layer < 3)
                        {
                            NTL::MulMod(prod[0], prod0[0], prod1[0], prod_bases[layer]);
                            if (k % 2 == 0)
                            {
                                if (layer == 2) 
                                    NTL::MulMod(prod[0], prod[0], p0p1a, prod_bases[layer]);
                                NTL::InvMod(prod[0], prod[0], prod_bases[layer]);
                            }
                        }

                        return prod;
                    },
                    layer == 3
                );
            prev_maps[k0].clear();
            prev_maps[k1].clear();

            std::vector<long> Pn {join_partitions(partitions, k*num_joined, (k+1)*num_joined)};
            next_partitions.push_back(std::move(Pn));
        }

        joined_partitions = std::move(next_partitions);
    }





}
