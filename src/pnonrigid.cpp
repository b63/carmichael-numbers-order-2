#include <iostream>
#include <vector>
#include <cstring>

#include <util.h>
#include <nonrigid.h>


int
main(int argc, char **argv)
{
    if (argc != 2)
        return 0;

    Factorization L { parse_factorization(argv[1]) };
    long max = 10000;

    //if (!parse_args(argc, argv, max, L))
    //    return 0;

    NTL::ZZ L_val;
    multiply_factors(L_val, L.primes, L.powers);
    std::cout << "L = " << L_val << "\n";

    std::vector<long> factors;
    get_nonrigid_primes(factors, L_val, max);
    printVec<long>(factors);
    std::cout << "\n";

    std::vector<std::array<long, 2> > pairs;
    get_gcd_Lfilter(pairs, factors, L_val);
    printVec<std::array<long, 2> >(pairs);



    std::cout << "\n";
}
