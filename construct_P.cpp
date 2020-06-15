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
#include <utility>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/RR.h>

#include "timer.h"
#include "util.h"
#include "counting_factors.h"


NTL::RR get_density(NTL::ZZ &L, size_t prime_size);
bool is_perfect_square(NTL::ZZ &n, const NTL::ZZ &n2);

void get_divisors(std::vector<NTL::ZZ> &divisors, 
        const std::vector<long> &prime_factors, const std::vector<long> &powers);

void get_divisors_with_factors(std::vector<NTL::ZZ> &divisors, 
        std::vector<std::vector<long>> &div_factors, std::vector<std::vector<long>> &div_powers,
        const std::vector<long> &prime_factors, const std::vector<long> &powers);

NTL::ZZ& multiply_factors(NTL::ZZ &prod, const std::vector<long> &factors, const std::vector<long> &powers);

NTL::ZZ construct_primes(
        std::map<const NTL::ZZ, std::vector<long> > &factor_map,
        const std::vector<long> &L_P_primes,
        const std::vector<long> &L_P_primes_powers,
        size_t SEIVE_SIZE = 1000
    );


NTL::ZZ construct_primes_2(
        std::map<const NTL::ZZ, std::vector<long> > &factor_map,
        const std::vector<long> &L_P_primes,
        const std::vector<long> &L_P_primes_powers,
        const long MAX
    );

NTL::ZZ populate_cofactor_map(std::map<const NTL::ZZ, std::vector<long> > &factor_map, 
        const std::vector<NTL::ZZ> &divisors, const std::vector<long> &primes);


void factor_n_1(std::vector<long>  &factors_n_1, std::vector<long> &powers_n_1,
        const std::vector<long> &factors_n2_1, const std::vector<long> &powers_n2_1,
        const NTL::ZZ &N);

void combine_factors( std::vector<long>  &factors, std::vector<long> &powers,
        const std::vector<long> &othr_factors, const std::vector<long> &othr_powers);




int main()
{
    // for timing
    init_timer();
    std::cout << std::fixed << std::setprecision(2);

    // naming note: L_P stands for L_PRIME -> L'
    // prime factorization of the L' parameter
    // prime factors  & corresponding powers
    const std::vector<long> L_P_PRIMES        {2, 5 }; 
    const std::vector<long> L_P_PRIMES_POWERS {2, 3 };

    // size of the sieve used to initially generate the prime numbers
    //const size_t SEIVE_SIZE = 10000000;

    // dictionary with co-factors k as key and vector of primes as values
    std::map<const NTL::ZZ, std::vector<long> > FACTOR_MAP;

    int time_id = start();
    construct_primes_2(FACTOR_MAP, L_P_PRIMES, L_P_PRIMES_POWERS, 1000);
    printTime(end(time_id));
}


NTL::ZZ construct_primes_2(
        std::map<const NTL::ZZ, std::vector<long> > &factor_map,
        const std::vector<long> &L_P_primes,
        const std::vector<long> &L_P_primes_powers,
        const long MAX
    ) 
{
    // for get_prime_factors
    init(MAX);

    // multiply the prime factors with appropriate powers to get L'
    NTL::ZZ L_P { 1 };
    multiply_factors(L_P, L_P_primes, L_P_primes_powers);
    std::cout << "L'=" << L_P << "\n";

    std::vector<NTL::ZZ> divisors;
    std::vector<std::vector<long>> factors;
    std::vector<std::vector<long>> powers;
    std::cout << "(generating divisors...)"<< std::flush;
    get_divisors_with_factors(divisors, factors, powers, L_P_primes, L_P_primes_powers);
    clrln(); std::cout << divisors.size() << " divisors\n";

    // skip the first and last divisor
    for (size_t i {1}; i+1 < divisors.size(); ++i)
    {
        NTL::ZZ &divisor { divisors[i] };
        NTL::ZZ multiple { divisor };
        NTL::ZZ N;

        for (long k {1}; multiple < MAX; ++k, multiple += divisor)
        {

            if (!is_perfect_square(N, multiple+1)) continue;
            const std::vector<long> &factors_k { get_prime_factors(k) };
            std::cout << "divisor " << divisor << ", N=" << N << " " ;
            printVec(factors_k);
            std::cout << "\n";

            std::vector<long> distinct_factors_k;
            std::vector<long> powers_k;
            collapse_factors(distinct_factors_k, powers_k, factors_k);
            combine_factors(factors[i], powers[i], distinct_factors_k, powers_k);

             printFactorization(factors[i], powers[i]);
            std::cout << "\n";
            std::vector<long> factors_n_1;
            std::vector<long> powers_n_1;
            factor_n_1(factors_n_1, powers_n_1, factors[i], powers[i], N);
            std::cout << "N^2-1 = " << N*N-1 << "\n"; 
            std::cout << "N-1 = ";
            printFactorization(factors_n_1, powers_n_1);
            std::cout << "\n\n";
            delete &factors_k;
        }
    }

    // temporary
    return L_P;
}

// does not assume the factor vector is sorted
void 
combine_factors( std::vector<long>  &factors, std::vector<long> &powers,
        const std::vector<long> &othr_factors, const std::vector<long> &othr_powers)
{
    const size_t len {factors.size()};
    const size_t othr_len {othr_factors.size()};

    for (size_t i {0}; i < othr_len; ++i)
    {
        bool common = false;
        long factor {othr_factors[i]};
        for (size_t j {0}; j < len; ++j)
        {
            if (factor == factors[j])
            {
                common = true;
                powers[j] += othr_powers[i];
                break;
            }
        }

        if (!common)
        {
            factors.push_back(factor);
            powers.push_back(othr_powers[i]);
        }
    }

}

void 
factor_n_1(std::vector<long>  &factors_n_1, std::vector<long> &powers_n_1,
        const std::vector<long> &factors_n2_1, const std::vector<long> &powers_n2_1,
        const NTL::ZZ &N)
{
    NTL::ZZ N_1 {N - 1};
    for(size_t i {0}; N_1 > 1 && i < factors_n2_1.size(); ++i)
    {
        long power {0};
        while (NTL::divide(N_1, N_1, factors_n2_1[i])) ++power;

        if (power > 0)
        {
            factors_n_1.push_back(factors_n2_1[i]);
            powers_n_1.push_back(power);
        }
    }
}



bool
is_perfect_square(NTL::ZZ &n, const NTL::ZZ &n2)
{
    NTL::ZZ root;
    NTL::SqrRoot(root, n2);

    NTL::ZZ square;
    NTL::sqr(square, root);

    while (square < n2)
    {
        root += 1;
        square = NTL::sqr(root);
    }

    if (square == n2)
    {
        n = root;
        return true;
    }
    else
    {
        return false;
    }
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
        std::map<const NTL::ZZ, std::vector<long> > &factor_map,
        const std::vector<long> &L_P_primes,
        const std::vector<long> &L_P_primes_powers,
        size_t sieve_size
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
    get_divisors(divisors, L_P_primes, L_P_primes_powers);
    clrln();

    std::cout << divisors.size() << " divisors\n";

    // get list of consecutive primes
    std::vector<long> primes;
    std::cout << "(generating primes...)" << std::flush;
    sieve_primes(primes, sieve_size);
    clrln();

    std::cout << "(filtering primes...)" << std::flush;
    // populate map with lists of primes
    const NTL::ZZ &k_max { populate_cofactor_map(factor_map, divisors, primes) };
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


/**
 * For every divisor in the vector `divisors`, the following process is repeated:
 *     for every p in the vector `primes` satisfying the modulo condition:
 *                  p^2 == 1 (mod divisor)
 *      the prime p is added to the map factor_map with co-factor k = (p^2-1)/divisor
 *      as the key
 *  the smallest co-factor k that accumulates the largest set of primes is returned
 */
NTL::ZZ
populate_cofactor_map(std::map<const NTL::ZZ, std::vector<long> > &factor_map, 
        const std::vector<NTL::ZZ> &divisors, const std::vector<long> &primes)
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
            std::vector<long> &modprimes {factor_map[k]};
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


/**
 * Populates vector `divisors` with divisors by going through every combination of prime factors in
 * `prime_factors` vector raised to every power up to corresponding entry in `powers` vector
 */
void get_divisors(std::vector<NTL::ZZ> &divisors, 
        const std::vector<long> &prime_factors, const std::vector<long> &powers)
{
    size_t primes = prime_factors.size();
    if (primes == 0)
        return;

    // stack of exponents for each prime factor
    int *stack {new int[primes]};
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

    delete[] stack;
}


/**
 * Populates vector `divisors` with divisors by going through every combination of prime factors in
 * `prime_factors` vector raised to every power up to corresponding entry in `powers` vector
 */
void get_divisors_with_factors(std::vector<NTL::ZZ> &divisors, 
        std::vector<std::vector<long>> &div_factors, std::vector<std::vector<long>> &div_powers,
        const std::vector<long> &prime_factors, const std::vector<long> &powers)
{
    size_t primes = prime_factors.size();
    if (primes == 0)
        return;

    // stack of exponents for each prime factor
    int *stack {new int[primes]};
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
                std::vector<long> f;
                std::vector<long> p;
                for (size_t i = 0; i < primes; ++i)
                {
                    int power = stack[i];
                    if (power > 0) 
                    {
                        prod *= std::pow(prime_factors[i], power);
                        f.push_back(prime_factors[i]);
                        p.push_back(power);
                    }
                }
                ++stack[top];
                // ignore 1
                divisors.push_back(prod);
                div_factors.push_back(f);
                div_powers.push_back(p);
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

    delete[] stack;
}


NTL::ZZ& 
multiply_factors(NTL::ZZ &prod, const std::vector<long> &factors, const std::vector<long> &powers)
{
    prod = 1;
    const size_t size { factors.size() };
    for (size_t i = 0; i < size; ++i)
    {
        prod *= std::pow(factors[i], powers[i]);
    }

    return prod;
}

