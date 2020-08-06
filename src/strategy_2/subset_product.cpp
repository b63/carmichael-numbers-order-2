#include <iostream>
#include <unordered_map>
#include <NTL/ZZ.h>

#include <util.h>
#include <config.h>
#include <strategy_2/subset_product.h>


typedef std::unordered_map<std::array<NTL::ZZ, 3>, std::vector<size_t>, ZZHasher<3>> map_type;

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

