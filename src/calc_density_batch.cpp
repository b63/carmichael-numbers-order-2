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

void
construct_primes(std::vector<long> &primes, const NTL::ZZ &L_val, const Factorization &L, long max = 0)
{
    NTL::ZZ p_max {NTL::SqrRoot(L_val+1)};
    std::cout << "upper bound on prime <= " << p_max << "\n";
    if (max)
        std::cout << "filtering primes <= " << max << " ...\n";
    else
        std::cout << "filtering primes <= " << p_max << " ...\n";

    NTL::PrimeSeq s;
    const size_t num_factors { L.primes.size() };
    size_t count = 0;
    size_t i { 0 };

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
                primes.push_back(p);
        }
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
        {
            primes.push_back(p);
        }
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
        std::cout << "reached small prime limit, starting from p=" << primes[primes.size()-1] << " ...\n";

        /* maxed out PrimeSeq, use NextPrime */
        NTL::ZZ p_zz {primes[primes.size()-1]};
        NextPrime(p_zz, p_zz);
        while (p_zz < p_max && (!max || p_zz <= max))
        {
            /* check to make sure p is not a factor of L */
            for(; i < num_factors && L.primes[i] < p_zz; ++i); /* NOP */
            if (i >= num_factors) 
                break;
            else if(NTL::divide(L_val, NTL::sqr(p_zz)-1))
                primes.push_back(NTL::conv<long>(p_zz));

            NextPrime(p_zz, p_zz);
#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
        }

        /* don't bother checking if p_zz is a factor of L */
        while (p_zz < p_max && (!max || p_zz <= max))
        {
            if(NTL::divide(L_val, NTL::sqr(p_zz)-1))
                primes.push_back(NTL::conv<long>(p_zz));

            NextPrime(p_zz, p_zz);

#if LOG_LEVEL >= 2
            if((count++ & STEP_MASK) == 0) std::cerr << "count: " << count << "\r";
#endif
        }

#if LOG_LEVEL >= 2
        std::cerr << "count: " << count << "\n";
#endif
    }
}


template <typename T>
void check_is_open(T stream, const char *name)
{
    if (!stream.is_open())
    {
        std::cerr << "Unable to open file \'" << name << "\'\n";
        return 0;
    }

    return 1;
}


int 
main(int argc, char **argv)
{
    const char *IN_FILE { argv[2] };
    const char *ALL_IN_FILE { argv[3] };
    const char *OUT_FILE { argv[4] };

    long max { NTL::conv<long>(argv[1]);
    std::ifstream in  {IN_FILE, std::ios::in};
    std::ifstream all_in  {ALL_IN_FILE, std::ios::in};
    std::ofstream out {OUT_FILE, std::ios::out | std::ios::trunc};

    if (!(check_is_open(in, IN_FILE) && check_is_open(all_in, ALL_IN_FILE)
            && check_is_open(out, OUT_FILE))
        return 1;


    std::string line;
    std::vector<long> primes;

    // first line contains the primes
    std::getline(in, line);
    size_t prev { 0 };
    while(1)
    {
        size_t pos = line.find(' ', prev);
        bool end = (pos == std::string::npos);
        std::string num { end ? line.substr(prev) : line.substr(prev, pos-prev) };
        if (num != " ")
            primes.push_back(NTL::conv<long>(num));

        if (end)
            break;
        else
            prev = pos+1;
    }

    while (std::getline(in, line))
    {
        size_t pos { line.find(':') };
        if (pos == std::string::npos) continue;
        std::string index_str { line.substr(0, pos) };
        std::string dist_str { line.substr(pos+1) };
    }
}


