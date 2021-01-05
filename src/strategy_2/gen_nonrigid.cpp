#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <NTL/ZZ.h>

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


Factorization
get_factorization(const NTL::ZZ &n)
{
    Factorization r;
    std::vector<NTL::ZZ> primes;
    factorize_slow(primes, r.powers, n);

    r.primes.reserve(primes.size());
    for (const auto &zz : primes) r.primes.push_back(NTL::conv<long>(zz));

    return r;
}


Factorization
combine(const Factorization &f1, const Factorization &f2)
{
    std::unordered_map<long,long> map;
    for (size_t i = 0, n = f1.primes.size(); i < n; i++)
        map[f1.primes[i]] += f1.powers[i];
    for (size_t i = 0, n = f2.primes.size(); i < n; i++)
        map[f2.primes[i]] += f2.powers[i];


    Factorization ff;
    for (auto it = map.cbegin(), cend = map.cend(); it != cend; ++it)
    {
        long prime = it->first, power = it->second;

        size_t i = 0, ff_size = ff.primes.size();
        for (; i < ff_size; ++i)
            if (prime < ff.primes[i])
                break;

        ff.primes.insert(ff.primes.cbegin()+i, prime);
        ff.powers.insert(ff.powers.cbegin()+i, power);
    }
    return ff;
}

/**
 * Format for commandline arguments(<> required, [] optional):
 *  [1]     <max> <limit> <p0> <p1> <min_terms> <max_terms> <L> <H>
 *  [2]     <max> <limit> <p0> <p1> <min_terms> <max_terms> <L> <H> < - | primes_path> [ -<min>,<max> | a_vals_path ]
 *
 * where
 *    max       - maximum prime when constructing set of primes P
 *    limit     - limit argument to pass to subroutine
 *    p0        - first non-rigid factor
 *    p1        - second non-rigid factor
 *    min_size  - minimum size subsets of P to consider
 *    max_size  - maximum size subsets of P to consider
 *    L         - parameter L specified as a prime factorization
 *               ex. "2^3 3^7" specifies L = 2^3 * 3^7
 *    H         - parameter H specified as an integer
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

    size_t argi = 1;
    const long max   { NTL::conv<long>(argv[argi++]) };
    const long limit { NTL::conv<long>(argv[argi++]) };
    const long p0    { NTL::conv<long>(argv[argi++]) };
    const long p1    { NTL::conv<long>(argv[argi++]) };
    const NTL::ZZ p0_zz {p0}, p1_zz {p1};
    const NTL::ZZ p02_1 {NTL::sqr(p0_zz)-1}, p12_1 {NTL::sqr(p1_zz)-1};

    const size_t min_terms { NTL::conv<size_t>(argv[argi++]) };
    const size_t max_terms { NTL::conv<size_t>(argv[argi++]) };

    Factorization L { parse_factorization(argv[argi++]) };
    Factorization p02_f { get_factorization(p02_1) };
    Factorization p12_f { get_factorization(p12_1) };
    Factorization LL_f {include_as_factor(L, include_as_factor(p02_f, p12_f))};

    const NTL::ZZ H {NTL::conv<NTL::ZZ>(argv[argi++])};

    /* print out information about parameter choices */
    NTL::ZZ L_val;
    multiply_factors(L_val, L.primes, L.powers);
    std::cout << "L = " << L_val << "\n";
    std::cout << "p0 = " << p0 << ", p1 = " << p1 << "\n";
    std::cout << "p0^2-1 = " << p02_1 << ", p1^2-1 = " << p12_1 << "\n";
    std::cout << "lcm(L, p0^2-1, p1^2-1) = " << LL_f << "\n";

    NTL::ZZ LL_val;
    multiply_factors(LL_val, LL_f.primes, LL_f.powers);
    std::cout << "                       = " << LL_val << "\n";


    /* read in primes from file or construct it */
    const std::array<long,2> nonrigid_factors { p0, p1 };
    std::vector<long> primes_set;
    if (argc >= 11 && strcmp(argv[argi], "-"))
    {
        read_primes_from_file(primes_set, argv[argi++]);
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
    if (argc >= 11 && argv[argi][0] != '-')
    {
        read_a_values_from_file(a_values, argv[argi++]);
    }
    else
    {
        char *comma = strchr(argv[argi], ',');
        if (comma == NULL) return 1;
        *comma = 0;
        int min {atoi(argv[argi++]+1)};
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

        gen_cprimes_4way_all(primes_set_a, nonrigid_factors, a_val, L_val, H, min_terms, max_terms, limit);
        std::cout << "\n\n";
    }
}
