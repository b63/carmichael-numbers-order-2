#include <iostream>
#include <vector>
#include <cstring>

#include <util.h>
#include <strategy_1/nonrigid.h>


int
main(int argc, char **argv)
{
    if (argc != 3)
        return 0;

    Factorization L { parse_factorization(argv[1]) };
    Factorization M { parse_factorization(argv[2]) };
    long max = 10000000;

    //if (!parse_args(argc, argv, max, L))
    //    return 0;

    NTL::ZZ L_val;
    multiply_factors(L_val, L.primes, L.powers);
    std::cout << "L = " << L_val << "\n";

    NTL::ZZ M_val;
    multiply_factors(M_val, M.primes, M.powers);
    std::cout << "M = " << M_val << "\n";

    std::vector<std::array<long, 2>> nonrigid_factors;
    get_nonrigid_factors(nonrigid_factors, L_val, M_val, max);
    printVec<std::array<long, 2> >(nonrigid_factors);
    std::cout << "\n";
}
