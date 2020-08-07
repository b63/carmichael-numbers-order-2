#include <iostream>
#include <unordered_map>
#include <NTL/ZZ.h>

#include <util.h>
#include <config.h>
#include <strategy_2/subset_product.h>


typedef std::unordered_map<std::array<NTL::ZZ, 3>, std::vector<size_t>, ZZHasher<3>> map_type;

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
    map_type map;

    /* calculate the total group size */
    const NTL::ZZ lcm {get_lcm<std::array<NTL::ZZ, num_bases> >(prod_base, num_bases)};

    /* split primes into two vectors */
    std::vector<long> h1_primes, h2_primes;
    split_half(h1_primes, h2_primes, primes);
    const size_t h1_primes_size {h1_primes.size()};

#if LOG_LEVEL >= 1
    std::cout << "|G| = " << lcm << "\n";
    std::cout << "storing inverse subset products of first half (size " << h1_primes_size << ") " << " ...\n";
    size_t inv_count {0};
#if LOG_LEVEL >= 2
    size_t count {0};
#endif
#endif
    /* go through every subset in first half, store inverse in hashmap */
    std::unique_ptr<std::vector<std::vector<size_t>>> subsets = subsetprod_mod<3>(h1_primes, prod_base, 
            [&](std::array<NTL::ZZ, num_bases>& prod_cache, 
                const std::vector<size_t> &indicies, size_t insert_index)->int
            {
                NTL::MulMod(prod_cache[2], prod_cache[2], p0p1a, prod_base[2]);
                std::array<NTL::ZZ, num_bases> invs;
                for(size_t i {0}; i < num_bases; ++i)
                   invs[i] = NTL::InvMod(prod_cache[i], prod_base[i]);
                std::vector<size_t> &vec {map[invs]};
                map[invs].push_back(insert_index);

#if LOG_LEVEL >= 1
                if(vec.size() == 1)
                    inv_count++;
#if LOG_LEVEL >= 2
                if((count++ & STEP_MASK) == 0)
                    std::cerr << std::setw(10) << "count: " << count << "\r";
#endif
#endif
                return 1;
            }, 1, 0);

#if LOG_LEVEL >= 1
#if LOG_LEVEL >= 2
        std::cerr << std::setw(10) << "count: " << count << "\r";
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
                const std::vector<size_t> ret {map[std::move(prod_cache)]};
                if (ret.size() > 0)
                {
                    for (auto ind : ret)
                    {
                        std::vector<long> cprimes;
                        std::vector<size_t> other_half {(*subsets)[ind]};
                        for(auto i : other_half) cprimes.push_back(h1_primes[i]);
                        for(auto i : indicies)   cprimes.push_back(h2_primes[i]);
                        std::cout << "found: " << cprimes << "\n";
                    }
                }
#if LOG_LEVEL >= 2
                if((count++ & STEP_MASK) == 0)
                    std::cerr << std::setw(10) << "count: " << count << "\r";
#endif
                return 0;
            }, 1, 0);

#if LOG_LEVEL >= 2
    std::cerr << std::setw(10) << "count: " << count << "\r";
#endif
}


void
gen_cprimes_2way_prob(
        const std::vector<long> &primes, 
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
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
    map_type map;

    /* calculate the total group size */
    const NTL::ZZ lcm {get_lcm<std::array<NTL::ZZ, num_bases> >(prod_base, num_bases)};

    /* split primes into two vectors */
    const size_t num_primes {primes.size()};
    std::vector<long> h1_primes, h2_primes;
    h1_primes.reserve(num_primes/2);
    for(auto it{primes.cbegin()},end=it+num_primes/2; it != end; ++it)
        h1_primes.push_back(*it);

    const size_t h1_primes_size {h1_primes.size()};
    const size_t h2_primes_size {num_primes-h1_primes_size};
    h2_primes.reserve(h2_primes_size);
    for(auto it{primes.cbegin()+h1_primes_size},end=primes.cend(); it != end; ++it)
        h2_primes.push_back(*it);

    const size_t sqrt_g {NTL::conv<size_t>(NTL::SqrRoot(lcm))};
    const size_t max_num_subsets {calc_max_subsets(h1_primes_size, 1, 0)};
    const size_t num_trials {MIN(sqrt_g, max_num_subsets<<1)};

#if LOG_LEVEL >= 2
    size_t count {0}; 
#if LOG_LEVEL >= 1
    std::cout << "|G| = " << lcm << "\n";
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
            std::cerr << std::setw(10) << "count: " << count << "\r";
#if LOG_LEVEL >= 1
        if (vec.size() == 1)
            inv_count++;
#endif
#endif
    }

#if LOG_LEVEL >= 2
    std::cerr << std::setw(10) << "count: " << subset_count << "\r";
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
        std::vector<size_t> &vec {map[std::move(prods)]};
        if (vec.size() > 0)
        {
            for (auto ind : vec)
            {
                std::vector<long> cprimes;
                std::vector<size_t> other_half {subsets[ind]};
                for(auto i : other_half) cprimes.push_back(h1_primes[i]);
                for(auto i : indicies)   cprimes.push_back(h2_primes[i]);
                std::cout << "found: " << cprimes << "\n";
            }
        }
#if LOG_LEVEL >= 2
        if((count++ & STEP_MASK) == 0)
            std::cerr << std::setw(10) << "count: " << count << "\r";
#endif
    }

#if LOG_LEVEL >= 2
    std::cerr << std::setw(10) << "count: " << count << "\r";
#endif
}
