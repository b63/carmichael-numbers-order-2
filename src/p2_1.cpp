#include <iostream>
#include <memory>
#include <vector>

#include <NTL/ZZ.h>
#include <util.h>
#include <counting_factors.h>


std::unique_ptr<std::vector<long>>
print_p2_1(size_t max)
{
    NTL::PrimeSeq seq;
    std::unique_ptr<std::vector<long>> all_primes_ptr {std::make_unique<std::vector<long>>()};

    long p {seq.next()};
    long max_prime {p};
    all_primes_ptr->push_back(p);
    while(p)
    {
        const size_t p2_1 {(size_t)(p*p-1)};
        if (p2_1 > max)
            break;

        std::vector<long> primes, powers;
        collapse_factors(primes, powers, *get_prime_factors(p*p-1));
        std::cout << p2_1 << " = ";
        for(size_t i{0}, len{primes.size()}; i < len; i++)
        {
            if (primes[i] > max_prime)
                max_prime = primes[i];

            if(i) std::cout << " ";
            std::cout << primes[i] << "^" << powers[i];
        }
        std::cout << "\n";

        p = seq.next();
        all_primes_ptr->push_back(p);
    }

    if (p && p < max_prime)
    {
        while (p && p < max_prime)
        {
            p = seq.next();
            all_primes_ptr->push_back(p);
        }
    }
    else if (p > max_prime)
    {
        std::unique_ptr<std::vector<long>> all_primes_cut_ptr {std::make_unique<std::vector<long>>()};
        for(auto prime : *all_primes_ptr)
        {
            if (prime > max_prime)
                break;
            all_primes_cut_ptr->push_back(prime);
        }

        return  all_primes_cut_ptr;
    }

    return all_primes_ptr;
}

int
main(int argc, char **argv)
{
    if (argc != 2) return 1;
    const size_t max {NTL::conv<size_t>(argv[1])};
    init(max);
    std::unique_ptr<std::vector<long>> all_primes_ptr {print_p2_1(max)};
    std::cout << "primes: ";
    for( auto p : *all_primes_ptr)
        std::cout << p << " ";
    std::cout << "\n";
}
