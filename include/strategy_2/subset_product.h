#ifndef S2_SUBSETPROD_H
#define S2_SUBSETPROD_H

#include <functional>
#include <vector>
#include <NTL/ZZ.h>

/* TEMPLATE FUNCTIONS */
template <unsigned long int N>
void
subsetprod_mod(const std::vector<long> &set, const std::array<NTL::ZZ, N> &bases,
        const std::function<void(std::array<NTL::ZZ, N>&, const std::vector<size_t>&)> &callback,
        size_t min_terms=2, size_t max_terms=0)
{
    constexpr size_t num_bases { N };

    const size_t set_size { set.size() };
    min_terms = MAX(min_terms, 2);
    max_terms = max_terms > 0 ? MIN(set_size, max_terms) : set_size;

    std::array<NTL::ZZ, num_bases> default_prod;
    std::fill(default_prod.begin(), default_prod.end(), 1);

    // go through all possible subset sizes starting from 2
    for (size_t factors {min_terms}; factors <= max_terms; factors++) 
    {
        std::vector<size_t> index_stack = {0};

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

                callback(mod_prod, index_stack);
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
}

#endif
