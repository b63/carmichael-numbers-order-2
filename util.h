#ifndef UTIL_H
#define UTIL_H
#include <vector>
#include <NTL/ZZ.h>

// clear `width` number of characters from terminal
// by writing that number of spaces followed by a
// carriage return
void clrln(size_t width=80);

// stores list of primes < size in `primes` vector
void sieve_primes(std::vector<long> &primes, size_t size);
// sieves `arr` , and `size` is lenght of arr
void sieve(long *arr, size_t size);


// TODO: finish implementation
bool poklington_test(NTL::ZZ N);

// uses NTL's built in prime generator to
// generate consecutive primes. Counts start from 0,
// so get_nth_prime(0) returns 2.
long get_nth_prime(size_t n);

#define UTIL_H_TEMPLATES 
#include "util.cpp"

#endif
