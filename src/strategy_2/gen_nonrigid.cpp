#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <fstream>
#include <algorithm>

#include <util.h>
#include <counting_factors.h>
#include <strategy_2/nonrigid.h>
#include <strategy_2/subset_product.h>


int
main(int argc, char **argv)
{
    if (argc <  5)
        return 0;

    const long max { NTL::conv<long>(argv[1]) };
    const long p0 { NTL::conv<long>(argv[2]) };
    const long p1 { NTL::conv<long>(argv[3]) };
    const long p02_1 { p0*p0 - 1};
    const long p12_1 { p1*p1 - 1};

    Factorization L { parse_factorization(argv[4]) };

    NTL::ZZ L_val;
    multiply_factors(L_val, L.primes, L.powers);
    std::cout << "L = " << L_val << "\n";

    std::cout << "p0 = " << p0 << ", p1 = " << p1 << "\n";
    std::cout << "p0^2-1 = " << p02_1 << ", p1^2-1 = " << p12_1 << "\n";

    const std::array<long,2> nonrigid_factors { p0, p1 };
    std::vector<long> primes_set;
    construct_primes_set(primes_set, nonrigid_factors, L_val, L, max);

    std::cout << "P = ";
    printVec<long>(primes_set);
    const size_t primes_set_size { primes_set.size() };
    std::cout << "\n|P| = " << primes_set_size << "\n";


    /* read-in/generate values for parameter a */
    std::vector<std::vector<long>> a_values;
    if (argc >= 6)
    {
        std::ifstream file { argv[5], std::ios::in};
        std::string line;
        while(std::getline(file, line))
        {
            std::vector<long> vals;
            if (parse_numbers(vals, line.c_str()) > 0)
                a_values.push_back(std::move(vals));
            else
                std::cout << "ignoring line, '" << line << "'\n";
        }
        file.close();
    }
    else
    {
        generate_a_values(a_values, primes_set, nonrigid_factors, 15, 0);
    }

    const size_t num_a { a_values.size() };
    std::cout << "number of a values: " << num_a << "\n\n";

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

        std::cout << "(|P| = " << primes_set_a.size() << ") trying a = " 
                << a_val << ", " << a_factors << "\n";
        gen_cprimes_2way_all(primes_set_a, nonrigid_factors, a_val, L_val, 1, 0);
        std::cout << "\n\n";
    }
}
