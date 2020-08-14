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
#include <memory>
#include <functional>


#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/RR.h>

#include <timer.h>
#include <util.h>
#include <counting_factors.h>
#include <construct_P.h>
#include <primality.h>
#include <config.h>

void test_method_2(const std::vector<long> &L_P_primes, const std::vector<long> &L_P_primes_powers, size_t lim)
{
    // multiply the prime factors with appropriate powers to get L'
    NTL::ZZ L_P { 1 };
    multiply_factors(L_P, L_P_primes, L_P_primes_powers);
    std::cout << "L'=" << L_P << "\nL'=";
    printFactorization(L_P_primes, L_P_primes_powers);
    std::cout << "\n";

    std::map<const long, std::vector<NTL::ZZ> > factor_map;

    int time_id = start();
    construct_primes_2(factor_map, L_P_primes, L_P_primes_powers, lim);
    printTime(end(time_id));

    // print stats
    std::cout << "\nstats:\n";
    print_stats<NTL::ZZ>(factor_map, L_P, 3);
}


void test_method_1(const std::vector<long> &L_P_primes, const std::vector<long> &L_P_primes_powers, size_t lim)
{
    // multiply the prime factors with appropriate powers to get L'
    NTL::ZZ L_P { 1 };
    multiply_factors(L_P, L_P_primes, L_P_primes_powers);
    std::cout << "L'=" << L_P << "\nL'=";
    printFactorization(L_P_primes, L_P_primes_powers);
    std::cout << "\n";

    std::map<const long, std::vector<long> > factor_map;

    int time_id = start();
    construct_primes(factor_map, L_P_primes, L_P_primes_powers, lim);
    printTime(end(time_id));

    // print stats
    std::cout << "\nstats:\n";
    print_stats<long>(factor_map, L_P, 3);
}


/********************** STANDARD METHOD ****************************/

/**
 * Simply goes through all primes p up to `max`, and checks if:
 *     (1) p does not divide L_val
 *     (2) p^2-1 divides L_val
 * If both conditions are satisfied, then added the primes to `primes` vector.
 * @param primes    vector<long> to which to add the primes that satify
 *                  conditions (1) and (2)
 * @param L_val     reference to NTL::ZZ value representing parameter L
 * @oaram L         Factorization object containing the prime factorization of L.
 *                  The vector of prime factors and powers are assumed to be sorted.
 *                  If they are not sorted, check for condition (1) might not work.
 */
void
construct_primes_standard(
        std::vector<long> &primes,
        const NTL::ZZ &L_val,
        const Factorization &L,
        size_t max)
{
    NTL::ZZ p_max {NTL::SqrRoot(L_val+1)};
    if (max && max < p_max)
        p_max = max;

#if LOG_LEVEL >= 1
    std::cout << "upper bound on prime <= " << NTL::SqrRoot(L_val+1) << "\n";
    std::cout << "filtering primes <= " << p_max << " ...\n";

#if LOG_LEVEL >= 2
    size_t count = 0;
#endif
#endif

    NTL::PrimeSeq s;
    const size_t num_factors { L.primes.size() };
    size_t i { 0 };

    long p = s.next();
    while ( p != 0 && p <= p_max)
    {
        /* check to make sure p is not a factor of L */
        for(; i < num_factors && L.primes[i] < p; ++i)
            ; /* NOP */

        if (i >= num_factors)
            break;
        else if (p != L.primes[i])
        {
            NTL::ZZ p2_1 { NTL::sqr(NTL::ZZ{p})-1 };
            if(NTL::divide(L_val, p2_1))
                primes.push_back(p);
        }
        p = s.next();
#if LOG_LEVEL >= 2
        if((count++ & STEP_MASK) == 0)
            std::cerr << "count: " << count << "\r";
#endif
    }

    /* NOTE: assming L does not contain a factor larger than what NTL::PrimeSeq supports */
    /* don't bother checking if p is a factor of L */
    while(p != 0 && p <= p_max)
    {
        if(NTL::divide(L_val, NTL::sqr(NTL::ZZ{p})-1))
        {
            primes.push_back(p);
        }
        p = s.next();
#if LOG_LEVEL >= 2
        if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
    }

#if LOG_LEVEL >= 2
    std::cerr << "count: " << count << "\n";
#endif

    if (p == 0 && p < p_max)
    {
        std::cout << "reached small prime limit, starting from p=" << primes[primes.size()-1] << " ...\n";

        /* maxed out PrimeSeq, use NextPrime */
        NTL::ZZ p_zz {primes[primes.size()-1]};
        NextPrime(p_zz, p_zz+1);
        while (p_zz < p_max)
        {
            /* check to make sure p is not a factor of L */
            for(; i < num_factors && L.primes[i] < p_zz; ++i); /* NOP */
            if (i >= num_factors) 
                break;
            else if(NTL::divide(L_val, NTL::sqr(p_zz)-1))
                primes.push_back(NTL::conv<long>(p_zz));

            NextPrime(p_zz, p_zz+1);
#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
        }

        /* don't bother checking if p_zz is a factor of L */
        while (p_zz < p_max)
        {
            if(NTL::divide(L_val, NTL::sqr(p_zz)-1))
                primes.push_back(NTL::conv<long>(p_zz));

            NextPrime(p_zz, p_zz+1);

#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) 
                std::cerr << "count: " << count 
                    << ", size: " << primes.size() << "\r";
#endif
        }

#if LOG_LEVEL >= 2
        std::cerr << "count: " << count << "\n";
#endif
    }
}


/********************** METHOD 2 ****************************/

/**
 * Constructs a set of primes using prime factorization 
 * of L' (L_P) parameter given through the `L_P_primes` and `L_P_primes_power`
 * lists.
 *
 * This method goes the multiples of each divisor of L' and
 * considers the ones that take th form:
 *                N^2 - 1  = divisor * k
 * Then runs the primality test of N to check if N is prime.
 * If N is prime then updates `factor_map` by adding N to the 
 * list of primes that the key k maps to
 *
 * ~~Returns the smallest co-factor k with the largest of list of primes.~~
 *
 * @param factor_map            reference to std::map that will be updated with 
 *                              co-factors as keys and list of corresponding primes
 *                              as values
 * @param L_P_primes            vector of prime factors of L'
 * @param L_P_primes_power      vector of powers of correspoding 
 *                              powers for the prime factors of L'
 * @param MAX                   maximum multiple of divisor to consider.
 *                              Only considers co-factors such that
 *                                          divisor*k <= MAX
 **/
void
construct_primes_2(
        std::map<const long, std::vector<NTL::ZZ> > &factor_map,
        const std::vector<long> &L_P_primes,
        const std::vector<long> &L_P_primes_powers,
        const long MAX
    )
{
    // for get_prime_factors
    init(MAX);

    std::vector<NTL::ZZ> divisors;
    std::vector<std::vector<long>> factors;
    std::vector<std::vector<long>> powers;

    get_divisors_with_factors(divisors, factors, powers, L_P_primes, L_P_primes_powers);
#if LOG_LEVEL >= 1
    std::cout << divisors.size() << " divisors\n";
#endif
    std::cout << "(filtering multiples of divisors...)\n";

    // skip the last divisor
    const size_t num_divisors {divisors.size()};
    for (size_t i {0}; i+1 < num_divisors; i++)
    {
        NTL::ZZ &divisor { divisors[i] };
        NTL::ZZ multiple { divisor };
        NTL::ZZ N;

        for (long k {1}; multiple < MAX; k++, multiple += divisor)
        {
            if (!is_perfect_square(N, multiple+1)) continue;
            const std::unique_ptr<std::vector<long>> factors_k { get_prime_factors(k) };

            std::vector<long> distinct_factors_k;
            std::vector<long> powers_k;
            collapse_factors(distinct_factors_k, powers_k, *factors_k);
            combine_factors(factors[i], powers[i], distinct_factors_k, powers_k);

            std::vector<long> factors_n_1;
            std::vector<long> powers_n_1;
            factor_n_1(factors_n_1, powers_n_1, factors[i], powers[i], N);

            bool prime { false };
            if (N == 2) /* gcd of N-1 will give error so seprate case for 2 */
            {
                prime = true;
            }
            else
            {
                /* test if N is composite */
                if (!quick_fermat_test(N)) continue;

                /* verify that N is prime */
                prime = pocklington(N, factors_n_1);
            }

            if (prime)
                factor_map[k].push_back(std::move(N));
        }
#if LOG_LEVEL >= 2
        if ((i&STEP_MASK) == 0)
            std::cerr << std::setw(10) << i << "/" << num_divisors << "\r" << std::flush;
#endif
    }
#if LOG_LEVEL >= 2
    std::cerr << std::setw(10) << num_divisors << "/" << num_divisors << "\n" << std::flush;
#endif
}


/**
 * Checks if Fermat's little theorem holds for a couple of small integers.
 * Returns false if a witness was found (i.e N is composite), true if no witness was found among the
 * integers checked (could be prime).
 * @param N  the integer to check if it is composite by trying to find a Fermat's witness for it
 */
bool
quick_fermat_test(const NTL::ZZ &N)
{
    static const std::vector<NTL::ZZ> witnesses {
        NTL::ZZ(2), NTL::ZZ(3), NTL::ZZ(4), NTL::ZZ( 5), 
        NTL::ZZ(6), NTL::ZZ(7), NTL::ZZ(8), NTL::ZZ(16)
    };
    const NTL::ZZ N_1 {N-1};

    NTL::ZZ mod;
    NTL::ZZ a;
    for(size_t i {0}; i < witnesses.size(); i++)
    {
        /* a and n must be relatively prime */
        NTL::GCD(mod, witnesses[i], N);
        if (!NTL::IsOne(mod)) continue;

        /* power mod gives error if any `a' is number if greater than base */
        NTL::rem(a, witnesses[i], N);

        NTL::PowerMod(mod, a, N_1, N);
        if (!NTL::IsOne(mod))
            return false;
    }

    return true;
}


/**
 * Populates vector `divisors` with divisors by going through every combination of prime factors in
 * `prime_factors` vector raised to every power up to corresponding entry in `powers` vector
 */
void
get_divisors_with_factors(std::vector<NTL::ZZ> &divisors, 
        std::vector<std::vector<long>> &div_factors, std::vector<std::vector<long>> &div_powers,
        const std::vector<long> &prime_factors, const std::vector<long> &powers)
{
#if LOG_LEVEL >= 2
    size_t count { 0 };
#endif
    size_t primes = prime_factors.size();
    if (primes == 0)
        return;

    // stack of exponents for each prime factor
    std::unique_ptr<int[]> stack = std::make_unique<int[]>(primes);
    size_t top { 0 };
    size_t max_top { primes - 1 };
    // keep track of the number of factors to call vector::reserve
    size_t last_size { 0 };

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
                f.reserve(last_size);
                p.reserve(last_size);

                for (size_t i = 0; i < primes; i++)
                {
                    int power = stack[i];
                    if (power > 0) 
                    {
                        prod *= std::pow(prime_factors[i], power);
                        f.push_back(prime_factors[i]);
                        p.push_back(power);
                    }
                }
                stack[top]++;
                // ignore 1
                if ((last_size = f.size()) > 0)
                {
                    divisors.push_back(std::move(prod));
                    div_factors.push_back(std::move(f));
                    div_powers.push_back(std::move(p));
#if LOG_LEVEL >= 2
                    if ((count++ & STEP_MASK) == 0) std::cerr << std::setw(10) << count << " divisors\r";
#endif
                }
            }
        }
        else
        {
            // expoenent for a prime factor is over max allowed
            if (top > 0) 
            {
                // pop current prime factor and increment exponent of the previous one
                stack[--top]++;
            } 
            else 
            {
                // reached bottom of stack
                break;
            }
        }
    }
#if LOG_LEVEL >= 2
    std::cerr << std::setw(10) << count << " divisors\n";
#endif
}


/**
 * Combines the factorization of two numbers (n1 and n2) into one but adding
 * factors to one of the vectors. Does not assume that the vectors are sorted.
 *
 * @param factors       reference to vector of prime factors of n1
 *                      prime factors of n2 are appended to this vector
 * @param powers        corresponding powers of the prime factors of n1,
 *                      will be modified as factors of n2 are appended to
 *                      `factors` vector
 * @param othr_factors  reference to vector of prime factors of n2   
 * @param othr_powers   reference to vector corresponding powers of
 *                      the prime factors of n2   
 */
void 
combine_factors( std::vector<long>  &factors, std::vector<long> &powers,
        const std::vector<long> &othr_factors, const std::vector<long> &othr_powers)
{
    const size_t len {factors.size()};
    const size_t othr_len {othr_factors.size()};

    for (size_t i {0}; i < othr_len; i++)
    {
        bool common = false;
        long factor {othr_factors[i]};
        for (size_t j {0}; j < len; j++)
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


/**
 * Extracts the factor of N-1 and stores it in `factors_n_1` and the powers
 * in `powers_n_1` given the prime factorization of N^2-1.
 *
 * @param factors_n_1   referece to vector where factors of N-1 should be stored
 * @param powers_n_1    referece to vector where corresponding powers of the factors
 *                      should be stored
 * @param factors_n2_1  reference to vector containing the prime factors of  N^2-1
 * @param powers_n2_1   referece to vector containing the corresponding powers of
 *                      the prime factors of N^2-1
 */
void 
factor_n_1(std::vector<long>  &factors_n_1, std::vector<long> &powers_n_1,
        const std::vector<long> &factors_n2_1, const std::vector<long> &powers_n2_1,
        const NTL::ZZ &N)
{
    NTL::ZZ N_1 {N - 1};
    for(size_t i {0}; N_1 > 1 && i < factors_n2_1.size(); i++)
    {
        long power {0};
        while (NTL::divide(N_1, N_1, factors_n2_1[i])) power++;

        if (power > 0)
        {
            factors_n_1.push_back(factors_n2_1[i]);
            powers_n_1.push_back(power);
        }
    }
}


/**
 * Returns true if `n2` is a perfect square, and stores the square root in `n`.
 * If `n2` is not a prefect square, the value of `n` is not changed.
 *
 * @param n  referce to a NTL::ZZ where the square root is stored if `n2` is a
 *           perfect square
 * @param n2 the number to check if it is a perfect square
 */
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





/********************** METHOD 1 ****************************/

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
 * ~~Returns the smallest co-factor k with the largest of list of primes.~~
 *
 * @param factor_map            reference to std::map that will be updated with 
 *                              co-factors as keys and list of corresponding primes
 *                              as values
 * @param L_P_primes            vector of prime factors of L'
 * @param L_P_primes_power      vector of powers of correspoding 
 *                              powers for the prime factors of L'
 * @param sieve_size            size of the sieve to use when generating the list
 *                              of primes so that only co-factors that satisfy
 *                                          divisor*k <= max_prime^2-1
 *                              are considered
 **/
void
construct_primes(
        std::map<const long, std::vector<long> > &factor_map,
        const std::vector<long> &L_P_primes,
        const std::vector<long> &L_P_primes_powers,
        size_t sieve_size
    )
{
    // get list of consecutive primes
    std::vector<long> primes;
    sieve_primes(primes, sieve_size);
#if LOG_LEVEL >= 1
    std::cout << "sieze size " << sieve_size << ", " << primes.size() << " primes\n";
#endif

    std::vector<NTL::ZZ> divisors;
    get_divisors(divisors, L_P_primes, L_P_primes_powers);
#if LOG_LEVEL >= 1
    std::cout << divisors.size() << " divisors\n";
    std::cout << "(filtering primes...)\n";
#endif

    // populate map with lists of primes
    populate_cofactor_map(factor_map, divisors, primes);
}


/**
 * For every divisor in the vector `divisors`, the following process is repeated:
 *     for every p in the vector `primes` satisfying the modulo condition:
 *                  p^2 == 1 (mod divisor)
 *      the prime p is added to the map factor_map with co-factor k = (p^2-1)/divisor
 *      as the key
 *  the smallest co-factor k that accumulates the largest set of primes is returned
 *
 *  NOTE: max{primes}^2 needs to fit in a long type
 */
void
populate_cofactor_map(
        std::map<const long, std::vector<long> > &factor_map, 
        const std::vector<NTL::ZZ> &divisors,
        const std::vector<long> &primes)
{
    const size_t num_primes   { primes.size() };
    const size_t num_divisors { divisors.size() };

    // skip last divisor
    for(size_t j = 0; j+1 < num_divisors; j++)
    {
        // TODO: maybe use gmp to do the modulo test
        const NTL::ZZ &divisor { divisors[j] };
        long ldivisor { NTL::conv<long>(divisor) };
        for (size_t i = 0;  i < num_primes; i++)
        {
            long prime = primes[i];

            NTL::ZZ_p::init(divisor);

            // calculate p^2 (mod divisor)
            NTL::ZZ_p p2 {prime};
            p2 *= prime;

            if (!NTL::IsOne(p2)) continue;

            // calculate co-factor k
            long k { (prime*prime -1)/ldivisor };

            //  add prime to map
            std::vector<long> &modprimes {factor_map[k]};
            modprimes.push_back(prime);
        }
#if LOG_LEVEL >= 2
        if ((j&STEP_MASK) == 0)
            std::cerr << std::setw(10) << j << "/" << num_divisors << "\r" << std::flush;
#endif
    }

#if LOG_LEVEL >= 2
    std::cerr << std::setw(10) << num_divisors << "/" << num_divisors << "\n";
#endif
}


/**
 * Populates vector `divisors` with divisors by going through every combination of prime factors in
 * `prime_factors` vector raised to every power up to corresponding entry in `powers` vector
 * Note: 1 is not included as a divisor
 */
void
get_divisors(std::vector<NTL::ZZ> &divisors, 
        const std::vector<long> &prime_factors, const std::vector<long> &powers)
{
#if LOG_LEVEL >= 2
    size_t count { 0 };
#endif
    size_t primes = prime_factors.size();
    if (primes == 0)
        return;

    // stack of exponents for each prime factor
    std::unique_ptr<int[]> stack = std::make_unique<int[]>(primes);
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
                for (size_t i = 0; i < primes; i++)
                {
                    int power = stack[i];
                    if (power > 0) 
                    {
                        prod *= std::pow(prime_factors[i], power);
                    }
                }
                stack[top]++;
                if (prod > 1)
                {
                    divisors.push_back(std::move(prod));
#if LOG_LEVEL >= 2
                    if((count++ & STEP_MASK) == 0)
                        std::cerr << std::setw(10) << count << " divisors...\r";
#endif
                }
            }
        }
        else
        {
            // expoenent for a prime factor is over max allowed
            if (top > 0) 
            {
                // pop current prime factor and increment exponent of the previous one
                stack[--top]++;
            } 
            else 
            {
                // reached bottom of stack
                break;
            }
        }
    }
#if LOG_LEVEL >= 2
    std::cerr << std::setw(10) << count << " divisors...\n";
#endif
}





