#ifndef COUNTING_FACTORS_H
#define COUNTING_FACTORS_H

#include <memory>

struct Product 
{
    long first_term { 0 };
    long second_term { 0 };
};

void init(size_t max);
std::unique_ptr<std::vector<long>> get_prime_factors(size_t n);
void collapse_factors(std::vector<long> &factors, std::vector<long> &powers, 
        const std::vector<long> &uncollapsed_factors);

#endif
