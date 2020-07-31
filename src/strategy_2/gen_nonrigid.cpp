#include <iostream>
#include <vector>
#include <cstring>
#include <algorithm>

#include <util.h>
#include <strategy_2/nonrigid.h>


int
main(int argc, char **argv)
{
    if (argc != 5)
        return 0;

    long max { NTL::conv<long>(argv[1]) };
    long p0 { NTL::conv<long>(argv[2]) };
    long p1 { NTL::conv<long>(argv[3]) };
    Factorization L { parse_factorization(argv[4]) };

    NTL::ZZ L_val;
    multiply_factors(L_val, L.primes, L.powers);
    std::cout << "L = " << L_val << "\n";

    std::cout << "p0 = " << p0 << ", p1 = " << p1 << "\n";

    const std::array<long,2> nonrigid_factors { p0, p1 };
    std::vector<long> primes_set;
    construct_primes_set(primes_set, nonrigid_factors, L_val, L, max);

    std::cout << "P = ";
    printVec<long>(primes_set);
    std::cout << "\nsize = " << primes_set.size() << "\n";

    std::vector<std::vector<long>> a_values;
    generate_a_values(a_values, primes_set, nonrigid_factors, 4);
    std::cout << "size: " << a_values.size() << "\n\n";


    const size_t num_a { a_values.size() };
    std::vector<std::vector<long>> cprimes;
    for(size_t i { 0 }; i < num_a; ++i)
    {
        NTL::ZZ a_val { 1 };
        std::vector<long> &a_factors { a_values[i] };
        for(auto it = a_factors.cbegin(), end=a_factors.cend();
                it != end; ++it)
            a_val *= *it;
        std::vector<long> primes_set_a;
        std::set_difference(primes_set.cbegin(), primes_set.cend(), 
                a_factors.cbegin(), a_factors.cend(),
                std::inserter(primes_set_a, primes_set_a.begin()));
        std::cout << primes_set_a.size() << ") trying a = " << a_val << ", " << a_factors << "\n";
        generate_cprimes(cprimes, primes_set, nonrigid_factors, a_val, a_factors, L_val, L);
    }
}
