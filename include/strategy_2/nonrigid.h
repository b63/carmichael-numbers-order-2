#ifndef S2_RIGID_H
#define S2_RIGID_H

#include <vector>
#include <util.h>


void generate_possible_factors(std::vector<std::array<long,2> > &factors, const NTL::ZZ &L_val, const long max=1000);

void generate_a_values(std::vector<std::vector<long> > &a_values, const std::vector<long> &primes,
        const std::array<long,2> &nonrigid_factors, size_t min_terms = 1, size_t max_terms=5);

void construct_primes_set(std::vector<long> &primes, const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &L_val, const Factorization &L, long max = 0);

void generate_cprimes(std::vector<std::vector<long>> &cprimes, const std::vector<long> primes, 
        const std::array<long, 2> &nonrigid_factors, 
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        size_t min_terms, size_t max_terms);
#endif
