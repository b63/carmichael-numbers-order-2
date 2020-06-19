#ifndef CONSTRUCT_P_H
#define CONSTRUCT_P_H


long construct_primes_2(
        std::map<const long, std::vector<NTL::ZZ> > &factor_map,
        const std::vector<long> &L_P_primes,
        const std::vector<long> &L_P_primes_powers,
        const long MAX
    );

void get_divisors_with_factors(std::vector<NTL::ZZ> &divisors, 
        std::vector<std::vector<long>> &div_factors, std::vector<std::vector<long>> &div_powers,
        const std::vector<long> &prime_factors, const std::vector<long> &powers);

void combine_factors( std::vector<long>  &factors, std::vector<long> &powers,
        const std::vector<long> &othr_factors, const std::vector<long> &othr_powers);

void factor_n_1(std::vector<long>  &factors_n_1, std::vector<long> &powers_n_1,
        const std::vector<long> &factors_n2_1, const std::vector<long> &powers_n2_1,
        const NTL::ZZ &N);

bool is_perfect_square(NTL::ZZ &n, const NTL::ZZ &n2);



long construct_primes(
        std::map<const long, std::vector<long> > &factor_map,
        const std::vector<long> &L_P_primes,
        const std::vector<long> &L_P_primes_powers,
        size_t SEIVE_SIZE = 1000
    );

long populate_cofactor_map(std::map<const long, std::vector<long> > &factor_map, 
        const std::vector<NTL::ZZ> &divisors, const std::vector<long> &primes);

void get_divisors(std::vector<NTL::ZZ> &divisors, 
        const std::vector<long> &prime_factors, const std::vector<long> &powers);



NTL::ZZ& multiply_factors(NTL::ZZ &prod, const std::vector<long> &factors, const std::vector<long> &powers);

NTL::RR get_density(NTL::ZZ &L, size_t prime_size);


#endif
