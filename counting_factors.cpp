#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "counting_factors.h"
#include "timer.h"
#include "util.h"

Product *MAP;
size_t PRODUCT_MAP_SIZE { 0 };


void init(size_t max)
{
    MAP =  new Product[max+1];
    PRODUCT_MAP_SIZE = max+1;

    size_t bound { static_cast<size_t>(ceil(sqrt(max+1))) };

    size_t p { 2 };
    while (p <= bound)
    {
        // skip if p is not prime
        if (MAP[p].first_term != 0) 
        {
            ++p;
            continue;
        }

        // keep track of second factor, sum = p * k
        size_t k { 1 };
        size_t sum { p };

        while (sum <= max)
        {
            Product &n = MAP[sum];
            // only write to MAP if it hasn't been written to before
            if (n.first_term == 0)
            {
                MAP[sum].first_term = (long) p;
                MAP[sum].second_term = (long) k;
            }

            sum += p;
            k += 1;
        }
    }
}

std::vector<long> &get_prime_factors(size_t n)
{
    if (n >= PRODUCT_MAP_SIZE)
    {
        throw std::invalid_argument(std::string("n too big, call init first")) ;
    }

    std::vector<long> *factors = new std::vector<long> {};

    while (true)
    {
        Product &prod = MAP[n];
        factors->push_back(prod.first_term);

        if (prod.second_term == 1)
        {
            // reached end of expansion
            break;
        }
        else
        {
            n = prod.second_term;
        }

    }

    return *factors;
}


/**
 * 'Collapses' the vector of factors given by `get_prime_factors`.
 * Separates the vector of factors including multiplicities into two vectors:
 * one containing the factors, and the other the powers.
 * The list of factors given must have repeated factors adjacent each other.
 * @param factors             reference to vector where distinct factors should be stored
 * @param powers              reference to vector where corresponding powers should be stored
 * @param uncollapsed_factors vector of factors that has repeated factors adjacent each other
 *                            eg. {2, 2, 4, 4, 5, 5, 5} but not {2, 3, 2, 3, 3, 4}
 */
void collapse_factors(std::vector<long> &factors, std::vector<long> &powers, 
        const std::vector<long> &uncollapsed_factors)
{
    const size_t len {uncollapsed_factors.size()};

    long prev = 0;
    long count = 0;
    for (size_t i {0}; i < len; ++i)
    {
        if (uncollapsed_factors[i] == prev)
        {
            ++count;
        }
        else 
        {
            if (prev != 0)
            {
                factors.push_back(prev);
                powers.push_back(count);
            }
            prev  = uncollapsed_factors[i];
            count = 1;
        }

    }
    // add the last factor to the vector, unless there were zero factors
    if (prev != 0)
    {
        factors.push_back(prev);
        powers.push_back(count);
    }
}
