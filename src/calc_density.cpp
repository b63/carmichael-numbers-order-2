#include <iomanip>
#include <iostream>
#include <vector>

#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include <timer.h>
#include <util.h>
#include <config.h>

void
construct_primes(std::vector<long> &primes, const NTL::ZZ &L_val, const Factorization &L, long max = 0)
{
    NTL::ZZ p_max {NTL::SqrRoot(L_val+1)};

#if LOG_LEVEL >= 1
    size_t count = 0;
    std::cout << "upper bound on prime <= " << p_max << "\n";
    if (max)
        std::cout << "filtering primes <= " << max << " ...\n";
    else
        std::cout << "filtering primes <= " << p_max << " ...\n";
#endif

    NTL::PrimeSeq s;
    const size_t num_factors { L.primes.size() };
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


int 
main(int argc, char **argv)
{
    // prime factors  & corresponding powers
    Factorization L;

    long factor, power;
    long max = 0;
    int limit { 0 };
    for (int i {1};  i < argc; i++)
    {
        if (*argv[i] == '-')
        {
            if (!limit && argv[i][1] == 'l')
            {
                size_t digits=0;
                char *ptr = argv[i]+2;
                for(; *ptr; digits++,ptr++)
                    if(*ptr < '0' || *ptr > '9') break;

                // at least 1 digits and no nondigit characters
                if (digits && !*ptr)
                {
                    max = NTL::conv<long>(argv[i]+2);
                    limit = 1;
                    continue;
                }
            }
        }
        parse_factor_string<long>(factor, power, argv[i]);
        L.primes.push_back(factor);
        L.powers.push_back(power);
    }

    if (L.primes.size() < 1)
    {
        return 0;
    }

    NTL::ZZ L_val;
    multiply_factors(L_val, L.primes, L.powers);
    std::cout << "L = " << L_val << "\n";

    std::vector<long> primes;
    init_timer();
    int id = start();
    construct_primes(primes, L_val, L, max);
    printTime(end(id));

    std::cout << "P = ";
    printVec<long>(primes);
    std::cout << "\n";

    NTL::RR density { get_density(L_val, primes.size()) };
    std::cout << "density = " << std::setprecision(10) << density << "\n";
}


