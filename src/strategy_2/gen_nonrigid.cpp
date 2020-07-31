#include <iostream>
#include <vector>
#include <cstring>

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
    generate_a_values(a_values, primes_set, nonrigid_factors, 10000);
    std::cout << "size: " << a_values.size() << "\n";
    printVec<std::vector<long>>(a_values);
}
