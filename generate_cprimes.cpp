#include <iostream>
#include <vector>
#include <cmath>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

void printVec(std::vector<int>&);
void seive(int*, int);
void printProd(std::vector<int>&, std::vector<int>&);

int main() {
	const int ORDER = 1;
	const int SEIVE_SIZE = 100000;
	// 27·33·52·7·11·13·17·19·29
	const int L_PRIMES[]        = {2, 3, 5,  7};
	const int L_PRIMES_POWERS[] = {3, 1, 1,  1};
	const int L_PRIMES_SIZE = sizeof(L_PRIMES)/sizeof(int);

	NTL::ZZ *L = new NTL::ZZ(1);
	for (int i = 0; i < L_PRIMES_SIZE; ++i)
	{
		*L *= std::pow(L_PRIMES[i], L_PRIMES_POWERS[i]);
	}
	NTL::ZZ_p::init(*L);
	std::cout << "L=" << *L << "\n";

	// create list of pimes
	int seive_array[SEIVE_SIZE];
	for (int i = 0; i < SEIVE_SIZE; ++i) {
		seive_array[i] = i+1;
	}

	seive(seive_array, SEIVE_SIZE);

	std::vector<int> primes;
	for (int i = 1; i < SEIVE_SIZE; ++i)
	{
		int n = seive_array[i];
		if (n != 0)
		{
			primes.push_back(n);
		}
	}

	std::cout << primes.size() << " primes\n\n";

	// filter list of primes
	std::vector<int> filtered_primes;
	int len = primes.size();
	for (int i = 0; i < len; ++i)
	{
		int prime = primes[i];
		bool include = true;
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
		
		filtered_primes.push_back(prime);
	}

	printVec(filtered_primes);
	std::cout << "\n" << filtered_primes.size() << " filtered primes" << "\n\n";

	// subset products of filetered_primes are possible carmicheal numbers
	// go over every subset of size factors
	int filtered_size = filtered_primes.size();
	std::vector<std::vector<int>> cprimes;

	for (int t=2; t <= filtered_size; ++t)
	{
		int factors = t;
		std::vector<int> index_stack = {0};
		std::cout << "\nchecking subsets of size " << t << " ... (found " << cprimes.size() << " till now)\n\n";

		while (index_stack.size() > 0)
		{
			int top_i = index_stack.size() - 1;
			// initializaiton
			if (index_stack[top_i] == -1)
			{
				index_stack[top_i] = index_stack.size() > 1 ? index_stack[top_i-1]+1 : 0;
			}

			if (index_stack[top_i] < filtered_size)
			{
				if (index_stack.size() < factors)
				{
					index_stack.push_back(-1);
					continue;
				}

				// check that this subset is a carmicheal prime
				bool carmicheal = true;
				NTL::ZZ prod(1);
				for (int j = 0; j < factors; ++j)
				{
					prod *= filtered_primes[index_stack[j]];
				}

				for (int j = 0; j < factors; ++j)
				{
					int p_factor = filtered_primes[index_stack[j]];
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

				if (carmicheal)
				{
					std::vector<int> copy(index_stack);
					cprimes.push_back(copy);
					printProd(index_stack, filtered_primes);
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
	
	int found = cprimes.size();
	std::cout << "found " << found << " carmicheal primes\n";
	for (int i = 0; i < found; ++i)
	{
		printProd(cprimes[i], filtered_primes);
		std::cout << "\n";
	}
}

void printProd(std::vector<int> &ind, std::vector<int> &arr)
{
	int factors = ind.size();
	for (int k = 0; k < factors; ++k)
	{
		if (k > 0) std::cout << " * ";
		std::cout << arr[ind[k]];
	}
}


void printVec(std::vector<int> &arr) 
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
