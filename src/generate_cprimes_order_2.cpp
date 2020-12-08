/**
 * Generates Carmicheal numbers of Order 2 using general Erdos construction.
 **/
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

#include <util.h>
#include <construct_P.h>
#include <subset_product.h>

/**
 * Format for commandline arguments:
 *       <max> <min_size> <max_size> <L>
 * where
 *    max      -  maximum prime when constructing set of primes P
 *    min_size - minimum size subsets of P to consider
 *    max_size - maximum size subsets of P to consider
 *    L        - parameter L specified as a prime factorization
 *               ex. "2^3 3^7" specifies L = 2^3 * 3^7
 */
int
main(int argc, char **argv)
{
    if (argc != 5)
        return 1;

    const size_t max { NTL::conv<size_t>(argv[1]) };
    size_t min_size { NTL::conv<size_t>(argv[2]) };
    size_t max_size { NTL::conv<size_t>(argv[3]) };
    Factorization L { parse_factorization(argv[4]) };

    NTL::ZZ L_val;
    multiply_factors(L_val, L.primes, L.powers);
    std::cout << "L = " << L_val << "\n";

    std::vector<long> primes_set;
    construct_primes_standard(primes_set, L_val, L, max);
    const size_t primes_set_size { primes_set.size() };
    std::cout << "P = " << primes_set << "\n";
    std::cout << "|P| = " << primes_set_size << "\n";

    const std::array<NTL::ZZ, 1> targets {NTL::ZZ{1}};
    const std::array<NTL::ZZ, 1> bases   { L_val };

    std::vector<std::vector<long>> cprimes;
    subsetprod_2way_all(cprimes, primes_set, targets, bases, min_size, max_size);

    for (auto v : cprimes)
    {
        std::cout << product<std::vector<long>, NTL::ZZ>(v) <<  " = " << v
            << "\n";
    }
}

