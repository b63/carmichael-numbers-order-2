#include <iostream>
#include <ostream>
#include <stdio.h>
#include <stdexcept>
#include <map>
#include <unordered_map>
#include <NTL/ZZ.h>

#include <util.h>

/*
 * Ensures that `n` includes a `factor` as a factor.
 * Note: Assumes that the vectors in `n` are sorted from
 * least to greatest.
 * @param n       reference to Factorization object consisting of two sorted vectors
 *                that contains n's prime factors and their powers
 * @param factor  factorization of the factor that `n` must have
 */
Factorization
include_as_factor(const Factorization &n, const Factorization &factor)
{
    const size_t factor_size = factor.primes.size();
    const size_t n_size      = n.primes.size();

    std::unordered_map<long, long> map {};

    size_t count = 0;
    for(size_t i = 0; i < n_size; i++)
    {
        long power = n.powers[i];
        if ( (map[n.primes[i]] += power) == power)
            count++; /* this prime hasn't been encountered before */
    }

    for(size_t i = 0; i < factor_size; i++)
    {
        long power = factor.powers[i];
        if ( (map[factor.primes[i]] += power) == power)
            count++; /* this prime hasn't been encountered before */
    }

    Factorization new_n;
    new_n.primes.reserve(count);
    new_n.powers.reserve(count);

    // TODO: test if fater than building unordered list, then sorting using
    //       algorithm of nlog(n) complexity?
    for (auto it = map.cbegin(), cend = map.cend(); it != cend; ++it)
    {
        long prime = it->first, power = it->second;

        size_t i = 0, nn_size = new_n.primes.size();
        for (; i < nn_size; ++i)
            if (prime < new_n.primes[i])
                break;

        new_n.primes.insert(new_n.primes.cbegin()+i, prime);
        new_n.powers.insert(new_n.powers.cbegin()+i, power);
    }

    // TODO: check for copy-elision
    return new_n;
}

std::ostream&
operator<<(std::ostream &os, const Factorization &f)
{
    const size_t len = f.primes.size();
    if (len != f.powers.size())
    {
        throw std::length_error("primes and powers vectors have unequal lengths");
    }

    for (size_t k = 0; k < len; k++)
    {
        if (k > 0) os << " * ";
        os << f.primes[k] << "^" << f.powers[k];
    }

    return os;
}


// generate list of primes using a sieve and store
// them in `primes` vector.
// @param primes vector instance where the primes will be stored
// @param size   size of the sieve used to generate primes
//               effectively limits the largest primes.
//               The largest prime will be <= size
void sieve_primes(std::vector<long> &primes, size_t size)
{
    // list of integers that will be sieved
    long *sieve_array = new long[size];
    for (size_t i = 0; i < size; i++) {
        sieve_array[i] = i+1;
    }

    sieve(sieve_array, size);

    // collect all the primes
    for (size_t i = 1; i < size; i++)
    {
        long n = sieve_array[i];
        if (n != 0)
        {
            primes.push_back(n);
        }
    }
}


// sieve an array of integers starting from 1 to create a  list of primes
// @prarm arr  pointer to the array, should contain [1, 2, ..., size-1]
// @param size the length of the array
void sieve(long *arr, size_t size)
{
    size_t index = 1;
    while (index < size)
    {
        long jump = arr[index];
        if (jump == 0)
        {
            index++;
            continue;
        }

        size_t  sieve_index = index + jump;

        for (; sieve_index < size; sieve_index += jump)
        {
            arr[sieve_index] = 0;
        }

        index++;
    }
}


bool poklington_test(const NTL::ZZ &N)
{
    static std::map<const NTL::ZZ, bool> _cache {};

    // for nonexistent key, _cache[key] is false
    // source: https://stackoverflow.com/questions/11058422/map-operator-and-bool-as-value
    if (_cache[N]) return true;

    const NTL::ZZ &N_1 {N-1};
    const NTL::ZZ &bound {NTL::SqrRoot(N)+1};

    // find f > sqrtN such that f | N-1
    NTL::ZZ f {1};
    NTL::ZZ r {N_1};
    std::vector<long> prime_factors {};
    for (size_t i {0}; f < bound; i++)
    {
        long ith_prime = get_nth_prime(i);
        // if ith_prime | r, then r /= ith_prime
        // otherwise r doesn't change
        if (NTL::divide(r, r, ith_prime))
        {
            prime_factors.push_back(ith_prime);
            f *= ith_prime;

            // keep dividing if possible
            while (f < bound && NTL::divide(r, r, ith_prime))
                f *= ith_prime;
        }
    }

    // TODO: need to ensure f and r are relatively prime
    return false;
}


/**
 * Returns the nth prime number using the PrimeSeq class
 * from the NTL library.
 */
long get_nth_prime(size_t n)
{
    static std::vector<long> _cache {0};
    static NTL::PrimeSeq s {};

    size_t size = _cache.size();
    if ( n < size)
    {
        return _cache[n];
    }

    _cache.reserve(n);
    while (size < n)
    {
        _cache.push_back(s.next());
        size++;
    }

    return _cache[n];
}


/** 
 * prints a carriage return ('\r') and `width` number
 * of spaces followed by another carriage return.
 * Helpful to 'erase' text written to a terminal
 */
void clrln(size_t width)
{
    char *spaces {new char[width+3]};
    spaces[0] = '\r';
    for (size_t i = 1; i <= width; i++) spaces[i] = ' ';
    spaces[width+1] = '\r';
    spaces[width+2] = '\0';
    std::cout << spaces;
    delete[] spaces;
}


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
    for (size_t i = 0; i < size; i++)
    {
        NTL::ZZ pow { NTL::power_ZZ(factors[i], powers[i]) };
        prod *= pow;
    }

    return prod;
}


/**
 * returns the density as a float, i.e 2^|prime_size| / L
 */
NTL::RR
get_density(NTL::ZZ &L, size_t prime_size)
{
    static const NTL::RR log2 { NTL::log(NTL::conv<NTL::RR>(2)) };
    NTL::RR mag_P { NTL::conv<NTL::RR>(prime_size) };
    NTL::RR l { NTL::conv<NTL::RR>(L) };
    if ( l <= 0)
        return NTL::conv<NTL::RR>(0);
    else
        return mag_P - NTL::log(L)/log2;
}

