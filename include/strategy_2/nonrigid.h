#ifndef S2_RIGID_H
#define S2_RIGID_H

#include <vector>
#include <functional>
#include <unordered_map>
#include <array>
#include <NTL/ZZ.h>

#include <util.h>
#include <subset_product.h>

template <size_t N>
using map_type = std::unordered_map<std::array<NTL::ZZ, N>, std::vector<std::vector<bool>>, ArrayHasher<NTL::ZZ, N>>;

void generate_possible_factors(std::vector<std::array<long,2> > &factors, const NTL::ZZ &L_val, const long max=1000);

void generate_a_values(std::vector<std::vector<long> > &a_values, const std::vector<long> &primes,
        const std::array<long,2> &nonrigid_factors, 
        const NTL::ZZ &L_val,
        size_t min_terms = 1, size_t max_terms=5);

void construct_primes_set(std::vector<long> &primes, const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &L_val, const Factorization &L, long max = 0);

NTL::ZZ&
calc_target_residue(NTL::ZZ &target, const NTL::ZZ &L_val,
        const NTL::ZZ &a_val,
        const std::array<long,2> &nonrigid_factors);

void gen_cprimes_2way_all(
        const std::vector<long> &primes, 
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        size_t min_size, size_t max_size);

void gen_cprimes_2way_prob(
        const std::vector<long> &primes, 
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        size_t num_trials,
        size_t min_size, size_t max_size);

map_type<1>*
gen_cprimes_4way_all(
        const std::array<std::vector<long>,4> &partition,
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        const NTL::ZZ &H, const NTL::ZZ &b,
        const std::function<std::array<NTL::ZZ,1>(const std::array<NTL::ZZ,1>&, const std::array<NTL::ZZ,1>&)> &callback,
        size_t min_size, size_t max_size,
        size_t LIMIT=0
    );

NTL::ZZ
get_group(
        const NTL::ZZ &L_val, const NTL::ZZ &a_val,
        const std::array<long,2> nonrigid_factors
    );

void
gen_cprimes_8way_pseudo(
        const std::array<std::vector<long>,4> &partition1,
        const std::array<std::vector<long>,4> &partition2,
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        const NTL::ZZ &H, const NTL::ZZ &sigma,
        size_t min_size, size_t max_size,
        size_t LIMIT=0
    );

void
gen_cprimes_8way_all(
        const std::vector<long> &primes,
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        size_t min_size, size_t max_size,
        size_t LIMIT=0
    );

void
gen_cprimes_8way_random(
        const std::vector<long> &primes,
        const std::array<long, 2> &nonrigid_factors,
        const NTL::ZZ &a_val, const NTL::ZZ &L_val,
        const std::array<NTL::ZZ, 4> &subgroups,
        size_t min_size, size_t max_size,
        size_t LIMIT=0
    );
#endif
