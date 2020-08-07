#ifndef S2_SUBSETPROD_H
#define S2_SUBSETPROD_H

#include <functional>
#include <vector>
#include <memory>
#include <NTL/ZZ.h>

#include <util.h>
#include <config.h>

void gen_cprimes_2way_all(
        const std::vector<long> &primes, 
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        size_t min_size, size_t max_size);

void
gen_cprimes_2way_prob(
        const std::vector<long> &primes, 
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        size_t min_size, size_t max_size);



/* TEMPLATE FUNCTIONS/CLASSES */

template <unsigned long int N>
struct ZZHasher
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

template <unsigned long int N>
std::unique_ptr<std::vector<std::vector<size_t>>>
subsetprod_mod(const std::vector<long> &set, const std::array<NTL::ZZ, N> &bases,
        const std::function<int(std::array<NTL::ZZ, N>&, const std::vector<size_t>&, size_t)> &callback,
        size_t min_terms=2, size_t max_terms=0)
{
    constexpr size_t num_bases { N };

    const size_t set_size { set.size() };
    min_terms = MAX(min_terms, 1);
    max_terms = max_terms > 0 ? MIN(set_size, max_terms) : set_size;

    std::array<NTL::ZZ, num_bases> default_prod;
    std::fill(default_prod.begin(), default_prod.end(), 1);

    // go through all possible subset sizes starting from 2
    std::unique_ptr<std::vector<std::vector<size_t>>> subsets_ptr 
        {std::make_unique<std::vector<std::vector<size_t>>>()};
    const size_t max_num_subsets {calc_max_subsets(set_size, min_terms, max_terms)};

#if LOG_LEVEL >= 1
    std::cout << "checking " << max_num_subsets << " subsets ...\n";
#endif

    size_t subset_count {0};
    subsets_ptr->reserve(max_num_subsets);
    for (size_t factors {min_terms}; factors <= max_terms; factors++) 
    {
        std::vector<size_t> index_stack {0};

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

            if (index_stack[top] < set_size) 
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

                int ret {callback(mod_prod, index_stack, subset_count)};
                if(ret == 1)
                {
                    subsets_ptr->push_back(index_stack);
                    subset_count++;
                }

                index_stack[top]++;
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

#if LOG_LEVEL >= 1
    std::cout << "stored  " << subsets_ptr->size() << " subsets \n";
#endif
    return subsets_ptr;
}

#endif
