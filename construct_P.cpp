/**
 * Explores the viability of construction a set of primes P from which to
 * create carmicheal numbers using the following algorithm:
 *	1. Pick a number L' (L_prime)
 *	2. For every divisor of L' find primes that satisfies d | p^2 - 1, and 
 *	   take note of the factor k = (p^2-1)/d
 *	3. Is there a factor k that comes up most often? If so let L = k*L'
 * Bobby
 **/
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

#include "timer.h"
#include "util.h"


void get_divisors(const std::vector<int> &prime_factors, const std::vector<int> &powers,
		std::vector<NTL::ZZ> &divisors);


// naming note: L_P stands for L_PRIME -> L'
// prime factorization of the L' parameter
const std::vector<int> L_P_PRIMES        {2, 3}; // prime factors
const std::vector<int> L_P_PRIMES_POWERS {2, 2}; // corresponding powers
const size_t L_P_PRIMES_SIZE = L_P_PRIMES.size();

// size of the seive used to initially generate the prime numbers
const size_t SEIVE_SIZE = 100;


int main() 
{
	// for timing
	init_timer();
	std::cout << std::fixed << std::setprecision(2);

	// multiply the prime factors with appropriate powers to get L'
	NTL::ZZ *L_P = new NTL::ZZ(1);
	for (size_t i = 0; i < L_P_PRIMES_SIZE; ++i)
	{
		*L_P *= std::pow(L_P_PRIMES[i], L_P_PRIMES_POWERS[i]);
	}
	 std::cout << "L=" << *L_P << "\n";

	std::vector<NTL::ZZ> divisors;
	get_divisors(L_P_PRIMES, L_P_PRIMES_POWERS, divisors);

	printVec(divisors);

	std::vector<int> primes;
	seive_primes(primes, SEIVE_SIZE);
	printVecInt(primes);
}




// populate divisors with divisors by going through every combination of prime factors in
// prime_factors vector raised to every power up to corresponding entry in powers vector
void get_divisors(const std::vector<int> &prime_factors, const std::vector<int> &powers, 
		std::vector<NTL::ZZ> &divisors)
{
	size_t primes = prime_factors.size();
	if (primes == 0)
		return;

	// stack of exponents for each prime factor
	int *stack = new int[primes];
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
				//std::cout << prod << "\n";
				divisors.push_back(prod);
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
}

