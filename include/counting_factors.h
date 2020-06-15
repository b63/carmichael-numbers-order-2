#ifndef COUNTING_FACTORS_H
#define COUNTING_FACTORS_H

struct Product 
{
    long first_term { 0 };
    long second_term { 0 };
};

void init(size_t max);
std::vector<long> &get_prime_factors(size_t n);

#endif
