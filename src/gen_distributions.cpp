#include <iostream>
#include <vector>
#include <cmath>


#include <util.h>

void gen_all_distributions(std::vector<std::vector<size_t> > &distributions, const std::vector<double> &weights, const double low, const double high);


bool
check_input(const char *ptr)
{
    size_t digits=0;
    for(; *ptr; digits++,ptr++)
        if(*ptr < '0' || *ptr > '9') break;

    // at least 1 digits and no nondigit characters
    if (digits && !*ptr)
        return true;
    else
        return false;
}

int
main(int argc, char **argv)
{
    constexpr static double log2 {log(2)};

    std::vector<double> weights;
    std::vector<long> primes;
    weights.reserve(argc);

    /* whether low/high was specified or not */
    int low = 0, high = 0;
    double dmax = 10.0, dmin = 9.0;
    char buf[50]; /* for formatting */

    for (int i {1};  i < argc; i++)
    {
        if (*argv[i] == '-')
        {
            if (!low && argv[i][1] == 'l')
            {
                // at least 1 digits and no nondigit characters
                if (check_input(argv[i]+2))
                {
                    dmin = NTL::conv<double>(argv[i]+2);
                    low = 1;
                    continue;
                }
            }
            else if (!high && argv[i][1] == 'h')
            {
                // at least 1 digits and no nondigit characters
                if (check_input(argv[i]+2))
                {
                    dmax = NTL::conv<double>(argv[i]+2);
                    high = 1;
                    continue;
                }
            }
        }


        if(check_input(argv[i]))
        {
            long p { NTL::conv<long>(argv[i]) };
            double logp { log(p)/log2 };
            weights.push_back(logp);
            primes.push_back(p);
        }
        else
        {
            snprintf(buf, 50, "invalid input, '%s'", argv[i]);
            throw std::invalid_argument(buf);
        }
    }

    if (weights.size() < 1)
        return 0;
    if (dmin > dmax)
    {
            snprintf(buf, 50, "low > high, %f>%f", dmin, dmax);
            throw std::invalid_argument(buf);
    }


    std::vector<std::vector<size_t> > distributions;
    gen_all_distributions(distributions, weights, dmin, dmax);

    for (size_t i = 0; i < primes.size(); i++)
    {
        if (i != 0) std::cout << ' ';
        std::cout << primes[i];
    }

    std::cout << "\n";
    size_t i {0};
    for (const auto &d : distributions)
    {
        std::cout << i << ": ";
        for(auto it=d.cbegin(); it!=d.cend(); it++ )
        {
            std::cout << *it;
            if (it+1 != d.cend())
                std::cout << " ";
        }
        std::cout << "\n";
        i++;
    }

}

/* TODO: do comments */
void
gen_all_distributions(std::vector<std::vector<size_t> > &distributions, 
        const std::vector<double> &weights, const double low, const double high)
{
    const size_t length { weights.size() };
    if (length < 1 || low >= high)
        return;

    // stack of exponents for each prime factor
    size_t *stack { new size_t[length] };
    double *sums { new double[length] };
    size_t top = 0;
    size_t max_top = length - 1;

    // go through every possible exponent for each prime factor
    while (true)
    {
        if (sums[top] < high)
        {

            if (top < max_top)
            {
                // every prime factor does not have exponent yet
                stack[top+1] = 0;
                sums[top+1] = sums[top];
                top++;
            }
            else
            {
                while (sums[top] < high)
                {
                    if (sums[top] >= low)
                    {
                        std::vector<size_t> powers;
                        powers.reserve(length);
                        for (size_t i = 0; i < length; i++)
                            powers.push_back(stack[max_top-i]);
                        distributions.push_back(std::move(powers));
                    }

                    sums[top]+= weights[max_top-top];
                    stack[top]++;
                }
            }
        }
        else
        {
            // expoenent for a prime factor is over max allowed
            if (top > 0) 
            {
                // pop current prime factor and increment exponent of the previous one
                top--;
                stack[top]++;
                sums[top] += weights[max_top-top];
            } 
            else 
            {
                // reached bottom of stack
                break;
            }
        }
    }

    delete[] stack;
    delete[] sums;
}
