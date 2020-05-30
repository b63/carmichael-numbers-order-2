#ifndef UTIL_H
#define UTIL_H
#include <vector>
#include <NTL/ZZ.h>

// prints the integers in arr specified by the array of indices ind 
// with a '*' as separator
void printProd(const std::vector<int> &ind, const std::vector<int> &arr);

// prints a vector of integers of type NTL:ZZ
void printVec(const std::vector<NTL::ZZ> &arr);

// prints a vector of integers
void printVecInt(const std::vector<int> &arr);

void seive(int *arr, size_t size);
void seive_primes(std::vector<int> &primes, size_t size);


#endif
