#include <iostream>
#include <stdexcept>
#include <map>

#include <NTL/ZZ.h>

#include "util.h"


void printProd(const std::vector<int> &ind, const std::vector<int> &arr)
{
    int factors = ind.size();
    for (int k = 0; k < factors; ++k)
    {
        if (k > 0) std::cout << " * ";
        std::cout << arr[ind[k]];
    }
}

void printVec(const std::vector<NTL::ZZ> &arr) 
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

void printVecInt(const std::vector<int> &arr) 
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


void seive_primes(std::vector<int> &primes, size_t size)
{
    // list of integers that will be seived
    int *seive_array = new int[size];
    for (size_t i = 0; i < size; ++i) {
        seive_array[i] = i+1;
    }
    
    seive(seive_array, size);
    
    // collect all the primes
    for (size_t i = 1; i < size; ++i)
    {
        int n = seive_array[i];
        if (n != 0)
        {
            primes.push_back(n);
        }
    }
}

// sieve an array of integers starting from 1 to create a  list of primes
// arr: [1, 2, 3, ..., size-1]
void seive(int *arr, size_t size)
{
    size_t index = 1;
    while (index < size)
    {
        int jump = arr[index];
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

long get_nth_prime(size_t n)
{
    static std::vector<long> _cache {0};
    static NTL::PrimeSeq s {};

    size_t size = _cache.size();
    if ( n < size)
    {
        return _cache[n];
    }

    while (size < n)
    {
        _cache.push_back(s.next());
        ++size;
    }

    return _cache[n];
}
