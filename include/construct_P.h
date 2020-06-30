#ifndef CONSTRUCT_P_H
#define CONSTRUCT_P_H


void construct_primes_2(
        std::map<const long, std::vector<NTL::ZZ> > &factor_map,
        const std::vector<long> &L_P_primes,
        const std::vector<long> &L_P_primes_powers,
        const long MAX
    );

bool quick_fermat_test(const NTL::ZZ &N);

void get_divisors_with_factors(std::vector<NTL::ZZ> &divisors, 
        std::vector<std::vector<long>> &div_factors, std::vector<std::vector<long>> &div_powers,
        const std::vector<long> &prime_factors, const std::vector<long> &powers);

void combine_factors( std::vector<long>  &factors, std::vector<long> &powers,
        const std::vector<long> &othr_factors, const std::vector<long> &othr_powers);

void factor_n_1(std::vector<long>  &factors_n_1, std::vector<long> &powers_n_1,
        const std::vector<long> &factors_n2_1, const std::vector<long> &powers_n2_1,
        const NTL::ZZ &N);

bool is_perfect_square(NTL::ZZ &n, const NTL::ZZ &n2);



void construct_primes(
        std::map<const long, std::vector<long> > &factor_map,
        const std::vector<long> &L_P_primes,
        const std::vector<long> &L_P_primes_powers,
        size_t seive_size = 1000
    );

void populate_cofactor_map(std::map<const long, std::vector<long> > &factor_map, 
        const std::vector<NTL::ZZ> &divisors, const std::vector<long> &primes);

void get_divisors(std::vector<NTL::ZZ> &divisors, 
        const std::vector<long> &prime_factors, const std::vector<long> &powers);


#endif
