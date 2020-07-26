#ifndef RIGID_H
#define RIGID_H

#include <vector>
#include <util.h>


void generate_nonrigid_cprimes(long p0, long p1, const Factorization &L_fact,
        const Factorization &M_fact, long sieve_size);
void subset_product_brute_force(std::vector<std::vector<size_t>> &cprimes, const std::vector<long> &primes);

#endif
