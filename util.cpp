#ifndef UTIL_H_TEMPLATES

#include <iostream>
#include <stdexcept>
#include <map>
#include <NTL/ZZ.h>

#include "util.h"


// generate list of primes using a sieve and store
// them in `primes` vector.
// @param primes vector instance where the primes will be stored
// @param size   size of the sieve used to generate primes
//               effectively limits the largest primes
void sieve_primes(std::vector<long> &primes, size_t size)
{
    // list of integers that will be sieved
    long *sieve_array = new long[size];
    for (size_t i = 0; i < size; ++i) {
        sieve_array[i] = i+1;
    }

    sieve(sieve_array, size);

    // collect all the primes
    for (size_t i = 1; i < size; ++i)
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
            ++index;
            continue;
        }

        size_t  sieve_index = index + jump;

        for (; sieve_index < size; sieve_index += jump)
        {
            arr[sieve_index] = 0;
        }

        ++index;
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
    for (size_t i {0}; f < bound; ++i)
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
        ++size;
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
    for (size_t i = 1; i <= width; ++i) spaces[i] = ' ';
    spaces[width+1] = '\r';
    spaces[width+2] = '\0';
    std::cout << spaces;
    delete[] spaces;
}


#else

// prints the integers in arr specified by the array of indices ind 
// with a '*' as separator
template <typename T>
void printProd(const std::vector<size_t> &ind, const std::vector<T> &arr)
{
    size_t factors = ind.size();
    for (size_t k = 0; k < factors; ++k)
    {
        if (k > 0) std::cout << " * ";
        std::cout << arr[ind[k]];
    }
}

// prints a vector of integers of type NTL:ZZ
template <typename T>
void printVec(const std::vector<T> &arr) 
{
    size_t size = arr.size();
    std::cout << "{";

    for (size_t i = 0; i < size; ++i)
    {
        if (i > 0) 
        {
            std::cout << ", ";
        }
        std::cout << arr[i];
    }

    std::cout << "}";
}


template <typename T>
void printFactorization(const std::vector<T> &factors, const std::vector<T> &powers)
{
    size_t len = factors .size();
    for (size_t k = 0; k < len; ++k)
    {
        if (k > 0) std::cout << " * ";
        std::cout << factors[k] << "^" << powers[k];
    }
}

#endif
