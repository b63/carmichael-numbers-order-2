#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include <timer.h>
#include <util.h>
#include <config.h>

#define PRE_LOAD 100

size_t
get_P_size(const NTL::ZZ &L_val, const Factorization &L, long max = 0)
{
    NTL::ZZ p_max {NTL::SqrRoot(L_val+1)};

#if LOG_LEVEL >= 1
    size_t count { 0 };
    std::cout << "upper bound on prime <= " << p_max << "\n";
    if (max)
        std::cout << "filtering primes <= " << max << " ...\n";
    else
        std::cout << "filtering primes <= " << p_max << " ...\n";
#endif

    NTL::PrimeSeq s;
    const size_t num_factors { L.primes.size() };
    size_t num_primes { 0 };
    size_t i { 0 };

    long p_prev { 0 };
    long p = s.next();
    while ( p != 0 && p <= p_max && (!max || p <= max))
    {
        /* check to make sure p is not a factor of L */
        for(; i < num_factors && L.primes[i] < p; ++i)
            ; /* NOP */

        if (i >= num_factors)
            break;
        else if (p != L.primes[i])
        {
            NTL::ZZ p2_1 { NTL::sqr(NTL::ZZ{p})-1 };
            if(NTL::divide(L_val, p2_1))
                num_primes++;
        }
        p_prev = p;
        p = s.next();
#if LOG_LEVEL >= 2
        if((count++ & STEP_MASK) == 0)
            std::cerr << "count: " << count << "\r";
#endif
    }

    /* NOTE: assming L does not contain a factor larger than what NTL::PrimeSeq supports */
    /* don't bother checking if p is a factor of L */
    while(p != 0 && p <= p_max && (!max || p <= max))
    {
        if(NTL::divide(L_val, NTL::sqr(NTL::ZZ{p})-1))
            num_primes++;
        p_prev = p;
        p = s.next();

#if LOG_LEVEL >= 2
        if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
    }

#if LOG_LEVEL >= 2
    std::cerr << "count: " << count << "\n";
#endif

    if (p == 0 && p < p_max)
    {
#if LOG_LEVEL >= 1
        std::cout << "reached small prime limit at p=" << p_prev << "\n";
#endif

        /* maxed out PrimeSeq, use NextPrime */
        NTL::ZZ p_zz {p_prev};
        NextPrime(p_zz, p_zz);
        while (p_zz < p_max && (!max || p_zz <= max))
        {
            /* check to make sure p is not a factor of L */
            for(; i < num_factors && L.primes[i] < p_zz; ++i); /* NOP */
            if (i >= num_factors) 
                break;
            else if(NTL::divide(L_val, NTL::sqr(p_zz)-1))
                num_primes++;

            NextPrime(p_zz, p_zz);
#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
        }

        /* don't bother checking if p_zz is a factor of L */
        while (p_zz < p_max && (!max || p_zz <= max))
        {
            if(NTL::divide(L_val, NTL::sqr(p_zz)-1))
                num_primes++;

            NextPrime(p_zz, p_zz);

#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
        }

#if LOG_LEVEL >= 2
        std::cerr << "count: " << count << "\n";
#endif
    }

    return num_primes;
}


template <typename T>
int check_is_open(const T &stream, const char *name)
{
    if (!stream.is_open())
    {
        std::cerr << "Unable to open file \'" << name << "\'\n";
        return 0;
    }

    return 1;
}


/**
 *
 */
int 
main(int argc, char **argv)
{
    if (argc != 5)
    {
        std::cerr << "not enough arguments\n";
        return 1;
    }

    const long max           { NTL::conv<long>(argv[1]) };
    const char *IN_FILE      { argv[2] };
    const char *ALL_IN_FILE  { argv[3] };
    const char *OUT_FILE     { argv[4] };

    std::ifstream in  {IN_FILE, std::ios::in};
    std::ifstream all_in  {ALL_IN_FILE, std::ios::in};
    std::ofstream out {OUT_FILE, std::ios::out | std::ios::trunc};

    if (!(check_is_open<std::ifstream>(in, IN_FILE) && check_is_open<std::ifstream>(all_in, ALL_IN_FILE)
            && check_is_open<std::ofstream>(out, OUT_FILE)))
    {
        return 1;
    }


    std::string line;
    std::vector<long> primes;

    // first line contains the primes
    std::getline(all_in, line);
    all_in.close();

    size_t prev { 0 };
    while(1)
    {
        size_t pos = line.find(' ', prev);
        bool end = (pos == std::string::npos);
        std::string num { end ? line.substr(prev) : line.substr(prev, pos-prev) };
        if (num != " ")
            primes.push_back(NTL::conv<long>(num.c_str()));

        if (end)
            break;
        else
            prev = pos+1;
    }

    const size_t num_primes { primes.size() };
    std::string lines_buffer[PRE_LOAD];

    while (in)
    {
        size_t i { 0 };
        for(; i < PRE_LOAD && std::getline(in, lines_buffer[i]); i++)
            ;

        for( size_t j { 0 }; j < i; ++j)
        {
            size_t pos { lines_buffer[j].find(':') };
            if (pos == std::string::npos)
                continue;

            std::string index_str { lines_buffer[j].substr(0, pos) };
            std::string dist_str { lines_buffer[j].substr(pos+1) };

            std::vector<long> powers;
            parse_numbers(powers, dist_str.c_str());

            if (num_primes != powers.size())
            {
                std::cerr << "unequal size, " << index_str << "\n";
                continue;
            }

            std::vector<long> primes_cpy { primes };
            Factorization L { std::move(primes_cpy), std::move(powers) };
            std::cout << index_str << ": " << L << "\n";

            NTL::ZZ L_val;
            multiply_factors(L_val, L.primes, L.powers);
            size_t P_size { get_P_size(L_val, L, max) };

            NTL::RR density { get_density(L_val, P_size) };
            
            out << index_str << ",";
            for(size_t k {0}; k < num_primes; ++k)
            {
                if (k) out << " ";
                out << L.powers[k];
            }
            out << dist_str << "," << density << "\n";
        }
    }

    in.close();
    out.close();
}


