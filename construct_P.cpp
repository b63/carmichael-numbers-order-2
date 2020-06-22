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
#include <functional>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/RR.h>

#include <timer.h>
#include <util.h>
#include <counting_factors.h>
#include <construct_P.h>

struct CoFactorSet
{
    /* number of primes associated with the set of cofactors */
    size_t num_primes;
    /* set of cofactors that have `num_primes` number of associated primes */
    std::vector<long> cofactors;
};

struct SizeCount
{
    /* number of primes */
    size_t num_primes;
    /* number of cofactors that have `num_primes` number of associated primes */
    size_t num_cofactors;
};


/**
 * Prints general statistics extracted from `factor_map`.
 * @param max_ranking    print the top `max_ranking` number of cofactors in terms of the
 *                       number of associated primes
*  @param T              type of the prime (long or NTL::ZZ)
 */
template <typename T>
void 
print_stats(const std::map<const long, std::vector<T>> &factor_map, size_t max_ranking)
{

    std::vector<CoFactorSet> ranking;
    std::vector<SizeCount> bag_sizes;

    for (auto it=factor_map.cbegin(); it != factor_map.cend(); ++it)
    {
        const long k      { it->first         };
        const size_t size { it->second.size() };

        /* update bag_sizes */
        size_t bag_size { bag_sizes.size() };
        size_t ind {0};
        bool found { 
            binary_search<std::vector<SizeCount>::const_iterator, SizeCount>
                (
                    ind, bag_sizes.begin(), 0, bag_size, 
                    [=](const SizeCount &item)->int 
                    {
                        if      (item.num_primes < size) return -1;
                        else if (item.num_primes > size) return  1;
                        else                             return  0;
                    }
                )
            };

        if (found)
        {
            bag_sizes[ind].num_cofactors += 1;
        }
        else
        {
            SizeCount entry {size, 1};
            if (ind == bag_size)
                bag_sizes.push_back(entry);
            else
                bag_sizes.insert(bag_sizes.begin()+ind, entry);
        }

        /* update ranking */
        found = binary_search<std::vector<CoFactorSet>::const_iterator, CoFactorSet>
            (
                ind, ranking.begin(), 0, ranking.size(),
                [=](const CoFactorSet &item)->int
                {
                    if      (item.num_primes < size) return -1;
                    else if (item.num_primes > size) return  1;
                    else                             return  0;
                }
            );
        if (found)
        {
            ranking[ind].cofactors.push_back(k);
        }
        else
        {
            CoFactorSet entry {size, std::vector<long>{k}};
            ranking.insert(ranking.begin()+ind, entry);
            if (ranking.size() >= max_ranking)
                ranking.erase(ranking.begin());
        }

    }

    /* print bag */
    size_t w {10};
    std::cout << std::setw(w) << "# primes" << "," << std::setw(w) << "# co-factors" << "\n";
    for(size_t i = 0; i < bag_sizes.size(); ++i)
    {
        const SizeCount &item {bag_sizes[bag_sizes.size() - i - 1]};
        std::cout << std::setw(w) << item.num_primes << ",";
        std::cout << std::setw(w) << item.num_cofactors << "\n";
    }
    std::cout << "\n";

    /* print ranking */
    if (ranking.size() > 0)
    {
        w = 5;
        for (auto it {ranking.crbegin()}; it != ranking.crend(); ++it)
        {
            std::cout << "num primes=" << it->num_primes << "\n";
            const std::vector<long> &cfs {it->cofactors};
            for (size_t j = 0; j < cfs.size(); ++j)
            {
                std::cout << "k=" << std::setw(w) << cfs[j] << " -> ";
                const std::vector<T> &primes {factor_map.at(cfs[j])};
                printVec<T>(primes);
                std::cout << "\n";
            }
            if (it+1 < ranking.crend())
                std::cout << "\n";
        }
    }
}


/**
 * Parses null-terminated `str` as "factor^power", and writes them to 
 * `factor` and `power`. If 'power' is not given, then 1 is assumed.
 * Throws std::invalid_argument if there were parsing errors.
 * Returns true if factor and power were successfully extracted.
 *
 * @param factor  reference to variable where factor is to be written
 * @param power   reference to variable where the power is to be written
 * @param str     pointer to first character of string to parse
 */
template <typename T>
bool
parse_factor_string(T &factor, T &power, char *str)
{
    size_t i { 0 };
    char *ch {str}, *mid {0};
    while(*ch)
    {
        if (*ch < '0' || *ch > '9')
        {
            /* character is not a digit */
            if (!i)
            {
                char buf[100];
                snprintf(buf, 100, "first character cannot be non-digit, '%s'", str);
                throw std::invalid_argument(buf);
            }
            else if (mid)
            {
                char buf[100];
                snprintf(buf, 100, "only digit characters allowed, '%s'", str);
                throw std::invalid_argument(buf);
            }
            else
            {
                mid = ch+1;
            }
        }
        ++i;
        ++ch;
    }
    if(!i) return false;

    if (!mid)
    {
        power = NTL::conv<T>(1);
        factor = NTL::conv<T>(str);
    }
    else
    {
        power = NTL::conv<T>(mid);
        char prev { *mid };
        *mid = 0;
        *mid = prev;
        factor = NTL::conv<T>(str);
    }

    return true;
}



void test_method_2(const std::vector<long> &L_P_primes, const std::vector<long> &L_P_primes_powers, size_t lim)
{
    // multiply the prime factors with appropriate powers to get L'
    NTL::ZZ L_P { 1 };
    multiply_factors(L_P, L_P_primes, L_P_primes_powers);
    std::cout << "L'=" << L_P << "\n";

    std::map<const long, std::vector<NTL::ZZ> > factor_map;

    int time_id = start();
    construct_primes_2(factor_map, L_P_primes, L_P_primes_powers, lim);
    printTime(end(time_id));

    // print stats
    std::cout << "\nstats:\n";
    print_stats<NTL::ZZ>(factor_map, 3);
}


void test_method_1(const std::vector<long> &L_P_primes, const std::vector<long> &L_P_primes_powers, size_t lim)
{
    // multiply the prime factors with appropriate powers to get L'
    NTL::ZZ L_P { 1 };
    multiply_factors(L_P, L_P_primes, L_P_primes_powers);
    std::cout << "L'=" << L_P << "\n";

    std::map<const long, std::vector<long> > factor_map;

    int time_id = start();
    construct_primes(factor_map, L_P_primes, L_P_primes_powers, lim);
    printTime(end(time_id));

    // print stats
    std::cout << "\nstats:\n";
    print_stats<long>(factor_map, 3);
}




int 
main(int argc, char **argv)
{
    // naming note: L_P stands for L_PRIME -> L'
    // prime factorization of the L' parameter
    // prime factors  & corresponding powers
    std::vector<long> L_P_primes  {}; 
    std::vector<long> L_P_primes_powers {};

    long factor, power;
    size_t max;

    int method {0}; // whether a method flag has been parsed or not
    int limit {0};  // whether a limit flag has been parsed or not
    for (int i {1};  i < argc; ++i)
    {
        if (*argv[i] == '-')
        {
            if (!method && argv[i][1] == 'm')
            {
                /* parse method option */
                char opt {argv[i][2]};
                if (opt == '1' || opt == '2')
                {
                    method = opt - '0';
                    continue;
                }
            }
            else if (!limit && argv[i][1] == 'l')
            {
                size_t digits=0;
                char *ptr = argv[i]+2;
                for(; *ptr; ++digits,++ptr)
                    if(*ptr < '0' || *ptr > '9') break;

                // at least 1 digits and no nondigit characters
                if (digits && !*ptr)
                {
                    max = NTL::conv<size_t>(argv[i]+2);
                    limit = 1;
                    continue;
                }
            }
        }


        parse_factor_string<long>(factor, power, argv[i]);
        L_P_primes.push_back(factor);
        L_P_primes_powers.push_back(power);
    }

    if (L_P_primes.size() < 1)
    {
        return 0;
    }

    printf("max %lu\n", max);

    // for timing
    init_timer();
    std::cout << std::fixed << std::setprecision(2);

    if (!method || method == 1)
        test_method_1(L_P_primes, L_P_primes_powers, max);
    else
        test_method_2(L_P_primes, L_P_primes_powers, max);
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

    //std::cout << "(generating divisors...)"<< std::flush;
    get_divisors_with_factors(divisors, factors, powers, L_P_primes, L_P_primes_powers);
    clrln(); std::cout << divisors.size() << " divisors\n";

    // skip the last divisor
    for (size_t i {0}; i+1 < divisors.size(); ++i)
    {
        NTL::ZZ &divisor { divisors[i] };
        NTL::ZZ multiple { divisor };
        NTL::ZZ N;

        for (long k {1}; multiple < MAX; ++k, multiple += divisor)
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

            factor_map[k].push_back(std::move(N));
        }
    }
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
    size_t primes = prime_factors.size();
    if (primes == 0)
        return;

    // stack of exponents for each prime factor
    int *stack {new int[primes]};
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
                if ((last_size = f.size()) > 0)
                {
                    divisors.push_back(std::move(prod));
                    div_factors.push_back(std::move(f));
                    div_powers.push_back(std::move(p));
                }
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
    std::cout << "sieze size " << sieve_size << ", " << primes.size() << " primes\n";

    std::vector<NTL::ZZ> divisors;
    get_divisors(divisors, L_P_primes, L_P_primes_powers);
    clrln(); std::cout << divisors.size() << " divisors\n";

    //clrln(); std::cout << "(filtering primes...)";
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
    for(size_t j = 0; j+1 < num_divisors; ++j)
    {
        // TODO: maybe use gmp to do the modulo test
        const NTL::ZZ &divisor { divisors[j] };
        long ldivisor { NTL::conv<long>(divisor) };
        for (size_t i = 0;  i < num_primes; ++i)
        {
            int prime = primes[i];

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
    }
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
                if (prod > 1)
                    divisors.push_back(std::move(prod));
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





/********************** UTIL ****************************/

/**
 * Multiplies each element of `factor` raised to the corresponding element of `powers`
 * and stores the final product in `prod`.
 *
 * @param prod    referece to NTL::ZZ where the final product is accumulated
 * @param factors vector of factors
 * @param powers  vector of powers for each factor in`factors`
 *
 * Returns a referce to `prod`.
 */
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

/**
 * returns the density as a float, i.e 2^|prime_size| / L
 */
NTL::RR
get_density(NTL::ZZ &L, size_t prime_size)
{
    // in arbitary precision, 2^|P| might be really big
    NTL::RR density {2};
    // density = density^primes_size
    NTL::pow(density, density, NTL::conv<NTL::RR>(prime_size));
    density /= NTL::conv<NTL::RR>(L);

    return density;
}

