#ifndef S1_RIGID_H
#define S1_RIGID_H

#include <vector>
#include <util.h>


void generate_nonrigid_cprimes(long p0, long p1, const Factorization &L_fact,
        const Factorization &M_fact, long sieve_size);
void subset_product_brute_force(std::vector<std::vector<size_t>> &cprimes, const std::vector<long> &primes);
void get_nonrigid_factors(std::vector<std::array<long, 2> > &nonrigid_factors, 
        const NTL::ZZ &L_val, const NTL::ZZ &M_val,
        size_t sieve_size);
void get_nonrigid_primes(std::vector<long> &primes, const NTL::ZZ &L, long max=10000);
void get_gcd_Lfilter(std::vector<std::array<long, 2> > &pairs, std::vector<long> &factors, const NTL::ZZ &L_val);


#endif
