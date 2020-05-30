#include <iostream>

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

