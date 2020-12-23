#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <fstream>
#include <algorithm>

#include <util.h>
#include <counting_factors.h>
#include <subset_product.h>
#include <strategy_2/nonrigid.h>


void
read_a_values_from_file(std::vector<std::vector<long> > &vec, const char *filename)
{
    std::ifstream file { filename, std::ios::in};
    std::string line;
    while(std::getline(file, line))
    {
        std::vector<long> vals;
        if (line.size() && line[0] != '#' && parse_numbers(vals, line.c_str()) > 0)
            vec.push_back(std::move(vals));
    }
    file.close();
}


void
read_primes_from_file(std::vector<long> &vec, const char *filename)
{
    std::ifstream file { filename, std::ios::in};
    std::string line;
    while(std::getline(file, line))
    {
        if (line.size() && line[0] != '#')
            vec.push_back(NTL::conv<long>(line.c_str()));
    }

    file.close();
}


/**
 * Format for commandline arguments(<> required, [] optional):
 *  [1]     <max> <p0> <p1> <min_terms> <max_terms> <L>
 *  [2]     <max> <p0> <p1> <min_terms> <max_terms> <L> < - | primes_path> [ -<min>,<max> | a_vals_path ]
 *
 * where
 *    max       - maximum prime when constructing set of primes P
 *    p0        - first non-rigid factor
 *    p1        - second non-rigid factor
 *    min_size  - minimum size subsets of P to consider
 *    max_size  - maximum size subsets of P to consider
 *    L         - parameter L specified as a prime factorization
 *               ex. "2^3 3^7" specifies L = 2^3 * 3^7
 *  For [2],
 *      < - | primes_path>   can be either '-' or a path to a file containing primes.
 *         If '-' or not specified, primes will be generated, otherwise primes constituting
 *         the primes set P will be read in from the path given.
 *
 *      [ - | a_vals_path]   can be either "-<min>,<max>" or path to file containing list of
 *         parameter a values to try. <min> and <max> refers to minimum subset to consideror
 *         where generating possible values for parameter a from  the primes set P.
 *         Since 2-set algorithm is used, subsets are drawn from halves.
 */
int
main(int argc, char **argv)
{
    if (argc <  7)
        return 0;

    const long max { NTL::conv<long>(argv[1]) };
    const long p0 { NTL::conv<long>(argv[2]) };
    const long p1 { NTL::conv<long>(argv[3]) };
    const NTL::ZZ p0_zz {p0}, p1_zz {p1};
    const NTL::ZZ p02_1 {NTL::sqr(p0_zz)-1}, p12_1 {NTL::sqr(p1_zz)-1};

    const size_t min_terms { NTL::conv<size_t>(argv[4]) };
    const size_t max_terms { NTL::conv<size_t>(argv[5]) };

    Factorization L { parse_factorization(argv[6]) };

    NTL::ZZ L_val;
    multiply_factors(L_val, L.primes, L.powers);
    std::cout << "L = " << L_val << "\n";

    std::cout << "p0 = " << p0 << ", p1 = " << p1 << "\n";
    std::cout << "p0^2-1 = " << p02_1 << ", p1^2-1 = " << p12_1 << "\n";

    const std::array<long,2> nonrigid_factors { p0, p1 };
    std::vector<long> primes_set;
    if (argc >= 9 && strcmp(argv[7], "-"))
    {
        read_primes_from_file(primes_set, argv[7]);
    }
    else
    {
        construct_primes_set(primes_set, nonrigid_factors, L_val, L, max);
        std::cout << "P = " << primes_set << "\n";
    }

    const size_t primes_set_size { primes_set.size() };
    std::cout << "|P| = " << primes_set_size << "\n";

    /* calculate and print total group size */
    const NTL::ZZ lcm {get_lcm<std::array<NTL::ZZ, 3> >(std::array<NTL::ZZ,3>{p02_1, p12_1, L_val}, 3)};
    const NTL::ZZ mod_G {eulers_toitent(lcm)};
    std::cout << "|G| = " << mod_G << " (" << ceil(NTL::log(mod_G)/log(2)) << " bits)\n";


    /* read-in/generate values for parameter a */
    std::vector<std::vector<long>> a_values;
    if (argc >= 9 && argv[8][0] != '-')
    {
        read_a_values_from_file(a_values, argv[8]);
    }
    else
    {
        char *comma = strchr(argv[8], ',');
        if (comma == NULL) return 1;
        *comma = 0;
        int min {atoi(argv[8]+1)};
        int max {atoi(comma+1)};
        *comma = ',';

        generate_a_values(a_values, primes_set, nonrigid_factors, L_val, min, max);
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
        gen_cprimes_4way_all(primes_set_a, nonrigid_factors, a_val, L_val, min_terms, max_terms);
        std::cout << "\n\n";
    }
}
