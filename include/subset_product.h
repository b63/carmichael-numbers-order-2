#ifndef S2_SUBSETPROD_H
#define S2_SUBSETPROD_H

#include <functional>
#include <vector>
#include <memory>
#include <cstdio>
#include <NTL/ZZ.h>

#include <util.h>
#include <config.h>

/* TEMPLATE FUNCTIONS/CLASSES */

template <typename T, size_t N>
struct ArrayHasher
{
    size_t
    operator()(const std::array<T, N> &arr) const
    {
        size_t hash {0};
        for (auto v : arr)
            hash ^= std::hash<T>{}(v)
                        + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        return hash;
    }
};


template <size_t N>
struct ArrayHasher<NTL::ZZ, N>
{
    size_t
    operator()(const std::array<NTL::ZZ, N> &arr) const
    {
        size_t hash {0};
        for (auto v : arr)
            hash ^= std::hash<unsigned long>{}(NTL::conv<unsigned long>(v))
                        + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        return hash;
    }
};

/**
 * Enumerate every single subset of `set` with sizes [min_terms, max_terms],
 * calculating the product of all the terms in the subset mod `bases`
 * (subset product) along the way. The subset-product, a copy of
 * the subset, and current subset count will be passed to the `callback` function.
 *
 * @param  set          vector of long integers from which to draw subsets
 * @param  bases        array of modulus bases to compute the subset-product
 * @param  callback     a function of that accepts the subset-product, the subset,
 *                      and the current subset count.
 * @param  min_terms    the minimum number of elements to draw
 * @param  max_terms    the maximum number of elements to draw
 * @param  LIMIT        maximum number of subsets to enumerate
 */
template <unsigned long int N>
void
subsetprod_mod(const std::vector<long> &set, const std::array<NTL::ZZ, N> &bases,
        const std::function<void(std::array<NTL::ZZ, N>&, const std::vector<size_t>&, size_t)> &callback,
        size_t min_terms=2, size_t max_terms=0, size_t LIMIT=0)
{
    constexpr size_t num_bases { N };

    const size_t set_size { set.size() };
    min_terms = MAX(min_terms, 1);
    max_terms = max_terms > 0 ? MIN(set_size, max_terms) : set_size;

    std::array<NTL::ZZ, num_bases> default_prod;
    std::fill(default_prod.begin(), default_prod.end(), 1);

    // go through all possible subset sizes starting from 2
    const size_t max_num_subsets {calc_max_subsets(set_size, min_terms, max_terms)};

#if LOG_LEVEL >= 1
    std::cout << "> enumerating " << max_num_subsets << " subsets [" << min_terms
        << "," << max_terms << "] (limit " << LIMIT << ") ...\n";
#endif

    size_t subset_count {0};
    for (size_t factors {min_terms}; factors <= max_terms; factors++)
    {
        std::vector<size_t> index_stack {0};
        index_stack.reserve(factors);

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

            if (index_stack[top] < set_size - (factors - index_size))
            {
                NTL::ZZ p_zz { set[index_stack[top]] };
                if (prod_size+1 < factors)
                {
                    std::array<NTL::ZZ, num_bases> prods;
                    for(size_t j { 0 }; j < num_bases; ++j)
                    {
                        NTL::ZZ prod { 1 };
                        if (prod_size)
                            prod = products[prod_top][j];
                        NTL::MulMod(prods[j], prod, p_zz, bases[j]);
                    }
                    products.push_back(std::move(prods));
                }

                if (index_size < factors)
                {
                    index_stack.push_back(index_stack[top]+1);
                    continue;
                }

                const std::array<NTL::ZZ, num_bases> &prods {
                    prod_size > 0 ? products[prod_top] : default_prod };

                /* multiply the current number p_zz to get the product
                 * of the current subset*/
                std::array<NTL::ZZ, num_bases> mod_prod;
                for(size_t j {0}; j < num_bases; ++j)
                    NTL::MulMod(mod_prod[j], prods[j], p_zz, bases[j]);

                callback(mod_prod, index_stack, subset_count);
                subset_count++;
                index_stack[top]++;

                if (LIMIT > 0 && subset_count >= LIMIT)
                {
                    return;
                }
            }
            else
            {
                index_stack.pop_back();
                if (top > 0) {
                    ++index_stack[top - 1];
                    products.pop_back();
                }

                continue;
            }
        }
    }
}


/**
 * Employ two-set subset-product algorithm to find subsets of `set` whose subset-product mod `bases`
 * is `targets`. Uses `subsetprod_mod` to go through all subsets of sizes [`min_size`, `max_size`]
 * from each halves of `set`.
 *
 * @param soln      vector of vector<long> to add solutions to
 * @param set       vector of long intergers to draw subsets from
 * @param targets   array of N residues (corresponds to the N values in `bases`)
 *                  that solutions should match
 * @param bases     array of N intergers that are the bases for the residues in `target`
 * @param min_size  minimum number of elements to draw from each half of `set`
 * @param max_size  maximum number of elements to draw from each half of `set`
 */
template <size_t N>
void subsetprod_2way_all(
    std::vector<std::vector<long>> &solns,
    const std::vector<long> &set,
    const std::array<NTL::ZZ, N>  &targets,
    const std::array<NTL::ZZ, N>  &bases,
    size_t min_size, size_t max_size)
{
    std::unordered_map<std::array<NTL::ZZ, N>, std::vector<std::vector<bool>>, ArrayHasher<NTL::ZZ, N>> map;

    /* split primes into two vectors */
    std::vector<long> first_half, second_half;
    split_half<long>(first_half, second_half, set);
    const size_t first_half_size {first_half.size()};

    /* calculate group size */
    const NTL::ZZ lcm {get_lcm<std::array<NTL::ZZ, N> >(bases, N)};
    const NTL::ZZ mod_G {eulers_toitent(lcm)};

    /* clip min_size and max_size */
    min_size = bound<size_t>(min_size, 1, first_half_size);
    if (max_size)
        max_size = bound<size_t>(max_size, min_size, first_half_size);
    else
        max_size = first_half_size;

    /* reserve space in hashmap */
    size_t num_subsets {calc_max_subsets(first_half_size, min_size, max_size)};
    constexpr size_t max_size_t {(size_t)-1};
    const size_t map_capacity {mod_G > max_size_t ? num_subsets
        : MIN(num_subsets, NTL::conv<size_t>(mod_G)) };

#if LOG_LEVEL >= 1
    std::cout << "|G| = " << mod_G << " (" << ceil(NTL::log(mod_G)/log(2)) << " bits)\n";
    std::cout << "subset sizes " << "[" << min_size << "," << max_size << "]: "
        << num_subsets << " subsets\n";
    std::cout << "map capacity: " << map_capacity << "\n";
    std::cout << "storing inverse subset products of first half (size "
        << first_half_size << ") " << " ...\n";
    size_t inv_count {0};
#if LOG_LEVEL >= 2
    size_t count {0};
#endif

    map.reserve(map_capacity);

#endif

    /* go through every subset in first half, store inverse in hashmap */
    subsetprod_mod<N>(first_half, bases,
            [&](std::array<NTL::ZZ, N>& prod_cache,
                const std::vector<size_t> &indicies, size_t subset_count)->void
            {
                std::array<NTL::ZZ, N> invs;
                for(size_t i{0}; i < N; ++i)
                {
                    if (!NTL::IsOne(targets[i]))
                        NTL::MulMod(prod_cache[i], prod_cache[i], targets[i], bases[i]);
                   invs[i] = NTL::InvMod(prod_cache[i], bases[i]);
                }

                std::vector<std::vector<bool>> &vec {map[invs]};

                // create vector<bool> representing this subset
                std::vector<bool> subset;
                subset.reserve(first_half_size);
                subset.resize(first_half_size, false);
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
            }, min_size, max_size);

#if LOG_LEVEL >= 1
#if LOG_LEVEL >= 2
    std::cerr << "count: " << count << "\r";
    count = 0;
#endif
    std::cout << "found " << inv_count << " distinct subset products\n";
    std::cout << "checking for inverse subset products of second half (size "
        << second_half.size() << ") ...\n";
#endif

    subsetprod_mod<N>(second_half, bases,
            [&](std::array<NTL::ZZ, N>& prod_cache,
                const std::vector<size_t> &indicies, size_t subset_count)->void
            {
                auto it { map.find(prod_cache)  };
                if (it != map.end())
                {
                    for (const std::vector<bool> &other_half : it->second)
                    {
                        std::vector<long> soln;
                        soln.reserve(indicies.size() + other_half.size());
                        for(size_t i {0}; i < other_half.size(); i++)
                            if (other_half[i])
                                soln.push_back(first_half[i]);

                        for(auto i : indicies)
                            soln.push_back(second_half[i]);
                        std::cout << "found: " << soln << "\n";
                        solns.push_back(std::move(soln));
                    }
                }
#if LOG_LEVEL >= 2
                if((count++ & STEP_MASK) == 0)
                    std::cerr << "count: " << count << "\r";
#endif
            }, min_size, max_size);
}


/**
 * Calculate the subset-product mod the M bases given in `bases`
 * for each of the subsets of `primes` given in`subsets` specified
 * as `vector<bool>`.
 *
 * @param dst       vector in which to store the subset-products
 * @param primes    vector<long> of elements from which subsets
 *                  in `subsets` are drawn
 * @param subsets   vector of subsets of `primes` specified as vector<bool>
 * @param baes      array of M bases for the subset-product
 */
template<unsigned long int M>
void
prod_mod_subsets(std::vector<std::array<NTL::ZZ, M>> &dst,
        const std::vector<long> &primes,
        const std::vector<std::vector<bool>> &subsets,
        const std::array<NTL::ZZ, M> &bases)
{
    const size_t size = subsets.size();
    size_t start = dst.size();
    dst.resize(start + size);

    for (size_t k = 0; k < size; k++)
    {
        std::array<NTL::ZZ, M> &prods = dst[start+k];
        const std::vector<bool> &subset = subsets[k];
        const size_t num_primes = subset.size();

        for (size_t i = 0; i < M; i++)
        {
            prods[i] = 1;

            for (size_t j = 0; j < num_primes; j++)
                if (subset[j])
                    NTL::MulMod(prods[i], prods[i], primes[j], bases[i]);
        }
    }

}


/**
 * Enumerate all subsets stored one of the maps (`src0` or `src1`) to find
 * keys (array of N subset-products) that are contained in both maps.
 *
 * If such pairs are found, then computes the new keys (array of M subset-products
 * in M bases) for the joined subset and stores them in `dst` using the key returned
 * by `callback.` The callback function is called with the products of each of the two
 * subsets to be joined in the new bases, and is expected to return the key that is used
 * to store the joined subset in map `dst`.
 *
 * The subsets are stored as `vector<bool>` so `primes0` is the parent set from which
 * the subsets stored in `src0` draw from; the same for `primes1` and the `vector<bool>`
 * subsets stored `src1`.
 *
 * @param dst           map to store the joined subsets using keys returned by `callback`
 * @param src0          one of the maps containing the subsets
 * @param src1          one of the maps containing the subsets
 * @param primes0       the parent set from which the subsets stored in `src0`
 *                      encoded as `vector<bool>` draw from
 * @param primes1       the parent set from which the subsets stored in `src1`
 *                      encoded as `vector<bool>` draw from
 * @param new_bases     array of M integers to use as the modulo base when
 *                      computing the new subset-product of the two subsets that are
 *                      passed on to `callback`
 * @param callback      function that returns the new key to use to store the joined subset
 *                      given the subset-product of the two subsets that are to be joined in
 *                      the new base
 * @param final_join    whether to not compute the new subset-product for each of the subsets
 *                      in the new base. If true, then this is the final so the new subset-product
 *                      are not computed and `callback` is not called to get the new key. The
 *                      default key of array of M elements is used as the key to store all
 *                      joined subsets.
 */
template <unsigned long int N, unsigned long int M>
void
join(std::unordered_map<std::array<NTL::ZZ, M>, std::vector<std::vector<bool>>, ArrayHasher<NTL::ZZ, M>> &dst,
        const std::unordered_map<std::array<NTL::ZZ, N>, std::vector<std::vector<bool>>, ArrayHasher<NTL::ZZ, N>> &src0,
        const std::unordered_map<std::array<NTL::ZZ, N>, std::vector<std::vector<bool>>, ArrayHasher<NTL::ZZ, N>> &src1,
        const std::vector<long> &primes0,
        const std::vector<long> &primes1,
        const std::array<NTL::ZZ,M> &new_bases,
        const std::function<std::array<NTL::ZZ,M>(const std::array<NTL::ZZ,M>&, const std::array<NTL::ZZ,M>&)> &callback,
        bool final_join=false
    )
{
#if LOG_LEVEL >= 1
    size_t pairs {0};
    size_t distinct_keys {0};
    size_t lookups {0};
    size_t count {0};
#endif
    /* default key, used if M == 0 or final_join == true */
    static const  std::array<NTL::ZZ, M> EMPTY_KEY;

    /* make map0 point to smaller of the two maps */
    auto *small_map = &src0;
    auto *big_map = &src1;
    bool flip {src1.size() < src0.size()};
    if (flip)
    {
        small_map = &src1;
        big_map   = &src0;
    }


    /* go through every item in map0 */
    const auto end = small_map->cend();
    for(auto it = small_map->cbegin(); it != end; it++)
    {
#if LOG_LEVEL >= 1
        if((lookups++ & STEP_MASK) == 0)
            fprintf(stderr, "\rlookups: %lu, pairs: %lu, subset count: %lu, keys: %lu",
                    lookups, pairs, count, distinct_keys);
#else
        lookups++;
#endif
        const std::array<NTL::ZZ, N> &key = it->first;

        auto it1 {big_map->find(key)};
        if (it1 != big_map->end())
        {
            /* get vector of subsets from src0 */
            const std::vector<std::vector<bool>> &subsets0 { flip ? src0.at(key) : it->second};
            /* get vector of subsets from src1 */
            const std::vector<std::vector<bool>> &subsets1 {!flip ? src1.at(key) : it->second};

            const size_t size0 = subsets0.size();
            const size_t size1 = subsets1.size();


            std::vector<std::array<NTL::ZZ, M>> products0, products1;
            if (M != 0 && !final_join)
            {
                /* store products of subsets in subsets0 and subset1
                * modulo new_bases in products0 and products1 */
                prod_mod_subsets(products0, primes0, subsets0, new_bases);
                prod_mod_subsets(products1, primes1, subsets1, new_bases);
            }
            else
            {
                /* final join, no need to compute the products of joined subset pairs
                 * in another base */
                std::cout << "key: ";
                for(auto &z : key) std::cout << z << " ";

                std::cout << "\n0: ";
                for (auto &s : subsets0)
                    print_bool_vec(s, primes0);

                std::cout << "\n1: ";
                for (auto &s : subsets1)
                    print_bool_vec(s, primes1);
                std::cout << "\n";
            }

            /* make a copy and store every combination of two subsets */
            for (size_t i = 0; i < size0; i++)
            {
                const std::vector<bool> &part0 {subsets0[i]};
                for (size_t j = 0; j < size1; j++)
                {
                    const std::vector<bool> &part1 {subsets1[j]};
                    std::vector<bool> joined {join_vector(part0, part1)};

                    std::vector<std::vector<bool>> *subsets_ptr {nullptr};

                    /* final join, use default key */
                    if (final_join || M == 0)  subsets_ptr = &dst[EMPTY_KEY];
                    /* get new key for joined subset using callback */
                    else  subsets_ptr = &dst[callback(products0[i], products1[j])];

                    subsets_ptr->push_back(std::move(joined));

#if LOG_LEVEL >= 1
                    if (subsets_ptr->size() == 1)  distinct_keys++;
#if LOG_LEVEL >= 2
                    if((count++ & STEP_MASK) == 0)
                        fprintf(stderr, "\rlookups: %lu, pairs: %lu, subset count: %lu, keys: %lu",
                                lookups, pairs, count, distinct_keys);
#else
                    count++;
#endif
#endif
                }
            }
#if LOG_LEVEL >= 1
            pairs++;
#endif
        }
    }

#if LOG_LEVEL >= 1
    fprintf(stderr, "\rlookups: %lu, pairs: %lu, subset count: %lu, keys: %lu\n",
            lookups, pairs, count, distinct_keys);
#endif
}


#endif

