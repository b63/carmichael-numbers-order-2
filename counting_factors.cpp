#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "counting_factors.h"
#include "timer.h"
#include "util.h"

Product *MAP;
size_t PRODUCT_MAP_SIZE { 0 };

int main()
{
    init_timer();
    int id { start() };

    size_t max { 1000000000 };

    start(id);
    init(max);
    printTime(end(id));

    long input;
    while ( true )
    {
        std::cout << "n = ";
        std::cin >> input;

        if (input < 2) break;

        start(id);
        const std::vector<long> factors = get_prime_factors(input);
        const time_metric &t = end(id);

        printVecLong(factors);
        std::cout << "\n";
        printTime(t);

        std::cout << std::endl;
    }
}

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
