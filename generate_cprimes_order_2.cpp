/**
 * Generates Carmicheal numbers of Order 2 using general Erdos construction.
 * Some optimizations specific to carmicheal numbers of order 2 may be utilized.
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

void filter_primes(const std::vector<long> &primes, std::vector<long> &dest);
void subset_product_brute_force(const std::vector<long> &,
                                std::vector<std::vector<size_t>> &);

// order of the carmicheal number
constexpr int ORDER = 2;
// size of the sieve used to initially generate the prime numbers
const size_t SEIVE_SIZE = 100000;

// prime factorization of the L parameter
const long L_PRIMES[]        = {2, 3, 5,  7}; // prime factors
const long L_PRIMES_POWERS[] = {3, 1, 1,  1}; // corresponding powers
const long L_PRIMES_SIZE = sizeof(L_PRIMES)/sizeof(long);


int main() {
    // for timing
    init_timer();
    std::cout << std::fixed << std::setprecision(2);

    // multiply the prime factors with appropriate powers to get L
    NTL::ZZ L { 1 };
    for (size_t i = 0; i < L_PRIMES_SIZE; ++i)
    {
        L *= std::pow(L_PRIMES[i], L_PRIMES_POWERS[i]);
    }
    NTL::ZZ_p::init(L);
    std::cout << "L=" << L << "\n";

    // list of integers that will be sieved
    long sieve_array[SEIVE_SIZE];
    for (size_t i = 0; i < SEIVE_SIZE; ++i) {
        sieve_array[i] = i+1;
    }

    sieve(sieve_array, SEIVE_SIZE);

    // collect all the primes
    std::vector<long> primes;
    for (size_t i = 1; i < SEIVE_SIZE; ++i)
    {
        long n = sieve_array[i];
        if (n != 0)
        {
            primes.push_back(n);
        }
    }

    std::cout << primes.size() << " primes\n";

    // filter list of primes and store in this vector
    std::vector<long> filtered_primes;
    filter_primes(primes, filtered_primes);

    std::cout << filtered_primes.size() << " filtered primes" << "\n";

    // subset products of filetered_primes are possible carmicheal numbers
    // go over every subset of size factors, store results in cprime list
    std::vector< std::vector<size_t> > cprimes;
    subset_product_brute_force(filtered_primes, cprimes);

    size_t found = cprimes.size();
    std::cout << "found " << found << " carmicheal primes\n";
    for (size_t i = 0; i < found; ++i)
    {
        printProd(cprimes[i], filtered_primes);
        std::cout << "\n";
    }
}


// filters the vector of primes using the following rule:
//     1. p must not divide L
//     2. p^2 - 1 must divide L
// the result is stored in the vector object dest
void filter_primes(const std::vector<long> &primes, std::vector<long> &dest)
{
    size_t len = primes.size();
    for (size_t i = 0; i < len; ++i)
    {
        long prime = primes[i];
        bool include = true;

        // check that p does not divide L
        for (size_t j = 0; j < L_PRIMES_SIZE; ++j)
        {
            if (L_PRIMES[j] == prime)
            {
                include = false;
                break;
            }
        }

        if (!include) continue;

        // check p^2 - 1 divides L
        NTL::ZZ_p p_mod_L(prime);
        p_mod_L *= prime;

        if (p_mod_L != 1) continue;
        
        dest.push_back(prime);
    }
}

// brute force approach to finding subset products with correct modulo.
// Searches all possible subsets of the vector primes and calculates the
// modulo product for each. 
//
// If the residue is what we want, the list of indicies
// is stored in the vector cprimes.
void subset_product_brute_force(const std::vector<long> &primes, std::vector<std::vector<size_t>> &cprimes)
{
    int timer_id = start(); // for timing with timer.h
    size_t num_primes = primes.size();
    // go through all possible subset sizes starting from 2
    for (size_t t = 2; t <= num_primes; ++t) 
    {
        size_t factors = t;
        std::vector<int> index_stack = {0};
        std::cout << "checking subsets of size " << t << " ...";
        std::cout << "(found " << cprimes.size() << " till now)\n";
        
        // start timer
        start(timer_id);
        // go through subsets of size factors suing a stack
        while (index_stack.size() > 0) 
        {
            size_t top_i = index_stack.size() - 1;
            // initializaiton
            if (index_stack[top_i] == -1) 
            {
                index_stack[top_i] = index_stack.size() > 1 ? index_stack[top_i - 1] + 1 : 0;
            }

            if ((size_t)index_stack[top_i] < num_primes) 
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

                // check that every prime factor satisfies the following divisibility
                // property
                //     p^2 -1 divides n
                for (size_t j = 0; j < factors; ++j) 
                {
                    long p_factor = primes[index_stack[j]];
                    NTL::ZZ p2(p_factor);
                    p2 *= p_factor;

                    if (prod % (p2 - 1) != 1) {
                    carmicheal = false;
                    break;
                    }
                }

                // create a copy of the indices for the prime numbers whose product
                // is a carmicheal number
                if (carmicheal) {
                    std::cout << "found ";
                    // TODO: use size_t for type fo index stack
                    std::vector<size_t> copy;
                    copy.reserve(index_stack.size());
                    for(size_t k=0; k < index_stack.size(); ++k) copy.push_back((size_t) index_stack[k]);
                    cprimes.push_back(copy);
                    printProd(copy, primes);
                    std::cout << "\n";
                }

                ++index_stack[top_i];
            } 
            else 
            {
                index_stack.pop_back();
                if (top_i > 0) {
                    ++index_stack[top_i - 1];
                }

                continue;
            }
        }

        // print elapsed time
        const time_metric &time { end(timer_id) };
        std::cout << "(cpu time " << time.cpu_time << " s)\n";
        std::cout << "(wall time " << time.wall_time << " s)\n";
    }
}


