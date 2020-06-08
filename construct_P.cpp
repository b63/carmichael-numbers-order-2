/**
 * Explores the viability of construction a set of primes P from which to
 * create carmicheal numbers using the following algorithm:
 *  1. Pick a number L' (L_prime)
 *  2. For every divisor of L' find primes that satisfies d | p^2 - 1, and 
 *     take note of the factor k = (p^2-1)/d
 *  3. Is there a factor k that comes up most often? If so let L = k*L'
 **/
#include <iomanip>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/RR.h>

#include "timer.h"
#include "util.h"


NTL::RR get_density(NTL::ZZ &L, size_t prime_size);
void get_divisors(const std::vector<int> &prime_factors, const std::vector<int> &powers,
        std::vector<NTL::ZZ> &divisors);

NTL::ZZ populate_cofactor_map(const std::vector<NTL::ZZ> &divisors, const std::vector<int> &primes,
        std::map<const NTL::ZZ, std::vector<int> > &factor_map);
NTL::ZZ construct_primes(
        const std::vector<int> &L_P_primes,
        const std::vector<int> &L_P_primes_powers,
        std::map<const NTL::ZZ, std::vector<int> > &factor_map,
        size_t SEIVE_SIZE = 1000
    );


void clrln()
{
    std::cout << "\r                                         \r";
}


int main()
{
    // for timing
    init_timer();
    std::cout << std::fixed << std::setprecision(2);

    // naming note: L_P stands for L_PRIME -> L'
    // prime factorization of the L' parameter
    // prime factors  & corresponding powers
    const std::vector<int> L_P_PRIMES        {2, 3, 5, 7, 11}; 
    const std::vector<int> L_P_PRIMES_POWERS {22, 13, 12, 4,  3};

    // size of the seive used to initially generate the prime numbers
    const size_t SEIVE_SIZE = 10000000;

    // dictionary with co-factors k as key and vector of primes as values
    std::map<const NTL::ZZ, std::vector<int> > FACTOR_MAP;

    int time_id = start();
    construct_primes(L_P_PRIMES, L_P_PRIMES_POWERS, FACTOR_MAP, SEIVE_SIZE);
    printTime(end(time_id));

}


/**
 * Constructs a set of primes using prime factorization 
 * of L' (L_P) parameter given through the `L_P_primes` and `L_P_primes_power`
 * lists..
 *
 * This method first constructs a set of primes using a sieve,
 * then goes through every divisor of L' and finds primes that
 * satisfy, divisor | (p^2-1). Each prime found if added to 
 * map `factor_map_ptr` with the co-factor k as the key.
 *
 * Returns the smallest co-factor k with the largest of list of primes.
 **/
NTL::ZZ construct_primes(
        const std::vector<int> &L_P_primes,
        const std::vector<int> &L_P_primes_powers,
        std::map<const NTL::ZZ, std::vector<int> > &factor_map,
        size_t seive_size
    ) 
{
    // multiply the prime factors with appropriate powers to get L'
    NTL::ZZ L_P { 1 };
    const size_t L_P_primes_size { L_P_primes.size() };
    for (size_t i = 0; i < L_P_primes_size; ++i)
    {
        L_P *= std::pow(L_P_primes[i], L_P_primes_powers[i]);
    }
    std::cout << "L'=" << L_P << "\n";

    std::vector<NTL::ZZ> divisors;
    std::cout << "(generating divisors...)"<< std::flush;
    get_divisors(L_P_primes, L_P_primes_powers, divisors);
    clrln();

    std::cout << divisors.size() << " divisors\n";

    // get list of consecutive primes
    std::vector<int> primes;
    std::cout << "(generating primes...)" << std::flush;
    seive_primes(primes, seive_size);
    clrln();

    std::cout << "(filtering primes...)" << std::flush;
    // populate map with lists of primes
    const NTL::ZZ &k_max { populate_cofactor_map(divisors, primes, factor_map) };
    const size_t k_max_size {  factor_map[k_max].size() };

    clrln();
    std::cout << "\nMaximum:\nk=" << k_max << ", size " << k_max_size << "\n";

    for (auto &item : factor_map)
    {
        const NTL::ZZ &k {item.first};
        size_t s {item.second.size()};
        if (s == k_max_size && k != k_max )
        {
            std::cout << "k=" << k << ", size " << s << "\n";
        }
    }

    // calculate density
    NTL::ZZ L {L_P * k_max};
    const NTL::RR &density = get_density(L, k_max_size);
    std::cout << "L=" << L << "\n";
    std::cout << "density " << std::setprecision(12) << density << "\n";

    return k_max;
}




// returns the density as a float, i.e 2^|prime_size| / L
NTL::RR get_density(NTL::ZZ &L, size_t prime_size)
{
    // in arbitary precision, 2^|P| might be really big
    NTL::RR density {2};
    // density = density^primes_size
    NTL::pow(density, density, NTL::conv<NTL::RR>(prime_size));
    density /= NTL::conv<NTL::RR>(L);

    return density;
}


// for every divisor in the vector `divisors`, the following process is repeated:
//     for every p in the vector `primes` satisfying the modulo condition:
//                  p^2 == 1 (mod divisor)
//      the prime p is added to the map factor_map with co-factor k = (p^2-1)/divisor
//      as the key
//  the smallest co-factor k that accumulates the largest set of primes is returned
NTL::ZZ populate_cofactor_map(const std::vector<NTL::ZZ> &divisors, const std::vector<int> &primes,
        std::map<const NTL::ZZ, std::vector<int> > &factor_map)
{
    const size_t num_primes   { primes.size() };
    const size_t num_divisors { divisors.size() };

    // to keep track of the co-factor k with max set of primes
    NTL::ZZ k_max;          // co-factor k
    size_t k_max_size {0};  // max size of set of primes

    // skip first and last divisor
    // they are 1 and L' itself
    for(size_t j = 1; j+1 < num_divisors; ++j)
    {
        const NTL::ZZ &divisor { divisors[j] };
        for (size_t i = 0;  i < num_primes; ++i)
        {
            int prime = primes[i];

            NTL::ZZ_p::init(divisor);

            // calculate p^2 (mod divisor)
            NTL::ZZ_p p2 {prime};
            p2 *= prime;

            if (!NTL::IsOne(p2)) continue;

            // calculate co-factor k
            NTL::ZZ k {prime};
            k *= prime;
            k -= 1;
            k /= divisor;

            //  add prime to map
            std::vector<int> &modprimes {factor_map[k]};
            modprimes.push_back(prime);

            // update max co-factor
            size_t modp_size = modprimes.size();
            if (modp_size > k_max_size)
            {
                k_max_size = modp_size;
                k_max = k;
            }
            else if (modp_size == k_max_size && k < k_max)
            {
                // if size of set of primes is the same, then pick
                // the smaller co-factor
                k_max = k;
            }
        }
    }

    return k_max;
}


// populate divisors with divisors by going through every combination of prime factors in
// prime_factors vector raised to every power up to corresponding entry in powers vector
void get_divisors(const std::vector<int> &prime_factors, const std::vector<int> &powers, 
        std::vector<NTL::ZZ> &divisors)
{
    size_t primes = prime_factors.size();
    if (primes == 0)
        return;

    // stack of exponents for each prime factor
    int *stack = new int[primes];
    size_t top = 0;
    size_t max_top = primes - 1;

    // go through every possible exponent for each prime factor
    stack[top] = 0;
    while (true)
    {
        if (stack[top] <= powers[top])
        {
            if (top < max_top)
            {
                // every prime factor does not have exponent yet
                stack[++top] = 0;
            }
            else
            {
                // every prime factor has a corresponding exponent
                // multiply it out to get the divisor
                NTL::ZZ prod{1};
                for (size_t i = 0; i < primes; ++i)
                {
                    int power = stack[i];
                    if (power > 0) 
                    {
                        prod *= std::pow(prime_factors[i], power);
                    }
                }
                ++stack[top];
                // ignore 1
                divisors.push_back(prod);
            }
        }
        else
        {
            // expoenent for a prime factor is over max allowed
            if (top > 0) 
            {
                // pop current prime factor and increment exponent of the previous one
                ++stack[--top];
            } 
            else 
            {
                // reached bottom of stack
                break;
            }
        }
    }
}

