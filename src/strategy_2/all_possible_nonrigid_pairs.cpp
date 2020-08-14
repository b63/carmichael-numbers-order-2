#include <iostream>
#include <vector>
#include <cstring>

#include <util.h>
#include <strategy_2/nonrigid.h>


/**
 * Prints out all possible pairs of primes that could be a viable
 * nonrigid factor for a given parameter choice L.
 * Namely, each factor in {p0, p1} must be units of the other
 * factor squared minus one.
 *                  gcd(p0, p1^2-1) == 1
 *                  gcd(p1, p0^2-1) == 1
 * Additionally, we must have p not dividing L and p^2-1 not dividing L.
 */
int
main(int argc, char **argv)
{
    if (argc != 3)
        return 0;

    long max = NTL::conv<long>(argv[1]);
    Factorization L { parse_factorization(argv[2]) };
    
    NTL::ZZ L_val;
    multiply_factors(L_val, L.primes, L.powers);
    std::cout << "L = " << L_val << "\n";

    std::vector<std::array<long, 2>> possible_pairs;
    generate_possible_factors(possible_pairs, L_val, max);
    printVec<std::array<long, 2> >(possible_pairs);
    std::cout << "\n";
}
