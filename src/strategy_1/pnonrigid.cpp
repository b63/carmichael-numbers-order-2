#include <iostream>
#include <vector>
#include <cstring>

#include <util.h>
#include <strategy_1/nonrigid.h>

/**
 * Takes the factorization of parameter L as the only cmdline argument,
 * and prints out the set of possible non-rigid factors up to a limit.
 * Then prints out the set of pairs of those factors that partly satisfy
 * the GCD contraints (under "partial: "). Then, finally the set of
 * pairs that satisfy all the GCD contraints imposed by L.
 */
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
