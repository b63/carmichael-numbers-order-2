#include <iomanip>
#include <iostream>
#include <vector>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/RR.h>

#include <timer.h>
#include <util.h>
#include <counting_factors.h>
#include <construct_P.h>
#include <primality.h>

int 
main(int argc, char **argv)
{
    // naming note: L_P stands for L_PRIME -> L'
    // prime factorization of the L' parameter
    // prime factors  & corresponding powers
    std::vector<long> L_P_primes  {}; 
    std::vector<long> L_P_primes_powers {};

    long factor, power;
    size_t max = 0;

    int method {0}; // whether a method flag has been parsed or not
    int limit {0};  // whether a limit flag has been parsed or not
    for (int i {1};  i < argc; i++)
    {
        if (*argv[i] == '-')
        {
            if (!method && argv[i][1] == 'm')
            {
                /* parse method option */
                char opt {argv[i][2]};
                if (opt == '1' || opt == '2')
                {
                    method = opt - '0';
                    continue;
                }
            }
            else if (!limit && argv[i][1] == 'l')
            {
                size_t digits=0;
                char *ptr = argv[i]+2;
                for(; *ptr; digits++,ptr++)
                    if(*ptr < '0' || *ptr > '9') break;

                // at least 1 digits and no nondigit characters
                if (digits && !*ptr)
                {
                    max = NTL::conv<size_t>(argv[i]+2);
                    limit = 1;
                    continue;
                }
            }
        }

        parse_factor_string<long>(factor, power, argv[i]);
        L_P_primes.push_back(factor);
        L_P_primes_powers.push_back(power);
    }

    if (L_P_primes.size() < 1)
    {
        return 0;
    }

    printf("max %lu\n", max);

    // for timing
    init_timer();
    std::cout << std::fixed << std::setprecision(2);

    if (!method || method == 1)
        test_method_1(L_P_primes, L_P_primes_powers, max);
    else
        test_method_2(L_P_primes, L_P_primes_powers, max*max);
}
