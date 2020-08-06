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

    /* split primes into two vectors */
    const size_t num_primes {primes.size()};
    std::vector<long> h1_primes, h2_primes;
    h1_primes.reserve(num_primes/2);
    for(auto it{primes.cbegin()},end=it+num_primes/2; it != end; ++it)
        h1_primes.push_back(*it);

    const size_t h1_primes_size {h1_primes.size()};
    h2_primes.reserve(num_primes-h1_primes_size);
    for(auto it{primes.cbegin()+h1_primes_size},end=primes.cend(); it != end; ++it)
        h2_primes.push_back(*it);

#if LOG_LEVEL >= 1
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
                std::array<NTL::ZZ, num_bases> invs;
                for(size_t i {0}; i < num_bases; ++i)
                   invs[i] = NTL::InvMod(prod_cache[i], prod_base[i]);
                const std::vector<size_t> ret {map[invs]};
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


size_t
calc_max_subsets(size_t set_size, size_t min_size, size_t max_size)
{
    size_t ret {0};
    max_size = max_size > 0 ? MIN(set_size, max_size) : set_size;
    for(; min_size <= max_size; min_size++)
        ret += binomial(set_size, min_size);
    return ret;
}

size_t binomial(unsigned long n, unsigned long m)
{
    if (m == 0 || n == m) return 1;
    else if (n == 0)      return 0;

    /* TODO: check the assembled output */
    /* allocate table  to store pascal's triangle */
    const size_t row {n+1};
    const size_t size {row*row};
    size_t *table {new size_t[size]};

    /* fill in pascal's triangle row by row */
    table[0] = 1;
    for(size_t i {1}; i <= n; ++i)
    {
        const size_t r {row*i};
        table[r] = i+1;
        table[r+i] = 1;
        for(size_t k{1}; k < i; ++k)
            table[r+k] = table[r-row+k] + table[r-row+k-1];
    }

    const size_t ret {table[row*(n-1)+m-1]};
    delete[] table;
    return ret;
}

