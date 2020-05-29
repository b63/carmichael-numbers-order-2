/**
 * Generates Carmicheal numbers using general Erdos construction.
 * Bobby
 **/
#include <iostream>
#include <vector>
#include <cmath>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

void printVec(const std::vector<int>&);
void printProd(const std::vector<int>&, const std::vector<int>&);

void seive(int *arr, int size);
void filter_primes(const std::vector<int> &primes, std::vector<int> &dest);
void subset_product(const std::vector<int>&, std::vector<std::vector<int> >&);


// order of the carmicheal number
const int ORDER = 1;
// size of the seive used to initially generate the prime numbers
const int SEIVE_SIZE = 100000;

// prime factorization of the L parameter
const int L_PRIMES[]        = {2, 3, 5,  7}; // prime factors
const int L_PRIMES_POWERS[] = {3, 1, 1,  1}; // corresponding powers
const int L_PRIMES_SIZE = sizeof(L_PRIMES)/sizeof(int);


int main() {
	// multiply the prime factors with appropriate powers to get L
	NTL::ZZ *L = new NTL::ZZ(1);
	for (int i = 0; i < L_PRIMES_SIZE; ++i)
	{
		*L *= std::pow(L_PRIMES[i], L_PRIMES_POWERS[i]);
	}
	NTL::ZZ_p::init(*L);
	std::cout << "L=" << *L << "\n";

	// list of integers that will be seived
	int seive_array[SEIVE_SIZE];
	for (int i = 0; i < SEIVE_SIZE; ++i) {
		seive_array[i] = i+1;
	}

	seive(seive_array, SEIVE_SIZE);

	// collect all the primes
	std::vector<int> primes;
	for (int i = 1; i < SEIVE_SIZE; ++i)
	{
		int n = seive_array[i];
		if (n != 0)
		{
			primes.push_back(n);
		}
	}

	std::cout << primes.size() << " primes\n";

	// filter list of primes and store in this vector
	std::vector<int> filtered_primes;
	filter_primes(primes, filtered_primes);

	std::cout << filtered_primes.size() << " filtered primes" << "\n\n";

	// subset products of filetered_primes are possible carmicheal numbers
	// go over every subset of size factors, store results in cprime list
	std::vector< std::vector<int> > cprimes;
	subset_product(filtered_primes, cprimes);

	int found = cprimes.size();
	std::cout << "found " << found << " carmicheal primes\n";
	for (int i = 0; i < found; ++i)
	{
		printProd(cprimes[i], filtered_primes);
		std::cout << "\n";
	}
}

// prints the integers in arr specified by the array of indices ind 
// with a '*' as separator
void printProd(const std::vector<int> &ind, const std::vector<int> &arr)
{
	int factors = ind.size();
	for (int k = 0; k < factors; ++k)
	{
		if (k > 0) std::cout << " * ";
		std::cout << arr[ind[k]];
	}
}


// prints a vector of integers
void printVec(const std::vector<int> &arr) 
{
	int size = arr.size();
	std::cout << "{";

	for (int i = 0; i < size; ++i)
	{
		if (i > 0) 
		{
			std::cout << ", ";
		}
		std::cout << arr[i];
	}

	std::cout << "}";
}

// sieve an array of integers starting from 1 to create a  list of primes
// arr: [1, 2, 3, ..., size-1]
void seive(int *arr, int size)
{
	int index = 1;
	while (index < size)
	{
		int jump = arr[index];
		if (jump == 0)
		{
			++index;
			continue;
		}

		int sieve_index = index + jump;

		for (; sieve_index < size; sieve_index += jump)
		{
			arr[sieve_index] = 0;
		}

		++index;
	}
}

// filters the vector of primes using the following rule:
//     1. p must not divide L
//     2. p^r - 1 must divide L for every r in 1,2,..,ORDER
// the result is stored in the vector object dest
void filter_primes(const std::vector<int> &primes, std::vector<int> &dest)
{
	int len = primes.size();
	for (int i = 0; i < len; ++i)
	{
		int prime = primes[i];
		bool include = true;

		// check that p does not divide L
		for (int j = 0; j < L_PRIMES_SIZE; ++j)
		{
			if (L_PRIMES[j] == prime)
			{
				include = false;
				break;
			}
		}

		if (!include) continue;

		// check p^r - 1 divides L
		NTL::ZZ_p p_mod_L(1);
		int r = 1;
		do { 
			p_mod_L *= prime;
			if (p_mod_L != 1)
			{
				include = false;
				break;
			}

			++r;
		} while (r <= ORDER);

		if (!include) continue;
		
		dest.push_back(prime);
	}
}

// brute force approach to finding subset products with correct modulo.
// Searches all possible subsets of the vector primes and calculates the
// modulo product for each. 
//
// If the residue is what we want, the list of indicies
// is stored in the vector cprimes.
void subset_product(const std::vector<int> &primes, std::vector< std::vector<int> > &cprimes)
{
	size_t num_primes = primes.size();
	// go through all possible subset sizes starting from 2
	for (size_t t=2; t <= num_primes; ++t)
	{
		unsigned int factors = t;
		std::vector<int> index_stack = {0};
		std::cout << "checking subsets of size " << t << " ...";
		std::cout <<  "(found " << cprimes.size() << " till now)\n";

		// go through subsets of size factors suing a stack
		while (index_stack.size() > 0)
		{
			size_t top_i = index_stack.size() - 1;
			// initializaiton
			if (index_stack[top_i] == -1)
			{
				index_stack[top_i] = index_stack.size() > 1 ? index_stack[top_i-1]+1 : 0;
			}

			if ((size_t) index_stack[top_i] < num_primes)
			{
				if (index_stack.size() < factors)
				{
					index_stack.push_back(-1);
					continue;
				}

				// check that this subset is a carmicheal prime
				bool carmicheal = true;
				NTL::ZZ prod(1);
				for (size_t j = 0; j < factors; ++j)
				{
					prod *= primes[index_stack[j]];
				}

				// check that every prime factor satisfies the following divisibility property
				//     p^r -1 divides n, where r ranges from 1,..,ORDER
				for (size_t j = 0; j < factors; ++j)
				{
					int p_factor = primes[index_stack[j]];
					NTL::ZZ p_r(1);
					for (int r = 1; r <= ORDER; ++r)
					{
						p_r *= p_factor;
						if (prod % (p_r-1) != 1)
						{
							carmicheal = false;
							break;
						}
					}
				}

				// create a copy of the indices for the prime numbers whose product
				// is a carmicheal number
				if (carmicheal)
				{
					std::vector<int> copy(index_stack);
					cprimes.push_back(copy);
					printProd(index_stack, primes);
					std::cout << "\n";
				}

				++index_stack[top_i];
			}
			else
			{
				index_stack.pop_back();
				if (top_i > 0)
				{
					++index_stack[top_i-1];
				}

				continue;
			}
		}
	}
}
