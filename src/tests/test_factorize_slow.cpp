#include <util.h>
#include <counting_factors.h>
#include <NTL/ZZ.h>


template <typename T1, typename T2>
bool
check_equal(std::vector<T1> v1, std::vector<T2> v2)
{
    const size_t s1 {v1.size()};
    const size_t s2 {v2.size()};
    if (s1 != s2) return false;

    for(size_t i = 0; i < s1; i++)
        if (v1[i] != v2[i])
            return false;
    return true;
}


int
test_example_1()
{
    long test_cases[] {2, 2<<8, 999, 123443};
    const size_t num_test_cases {sizeof(test_cases)/sizeof(test_cases[0])};

    long max {1};
    for (size_t i {0}; i < num_test_cases; i++)
        if (test_cases[i] > max) max = test_cases[i];
    init(max);

    int ok {1};
    for (size_t i {0}; i < num_test_cases; i++)
    {
        std::vector<long> primes, powers;
        collapse_factors(primes, powers, *get_prime_factors(test_cases[i]));

        std::vector<NTL::ZZ> primes_test;
        std::vector<long> powers_test;
        factorize_slow(primes_test, powers_test, NTL::ZZ{test_cases[i]});

        if (!check_equal<NTL::ZZ, long>(primes_test, primes)
            || !check_equal<long, long>(powers_test, powers))
        {
            ok = 0;
            std::cerr << "FAIL case " << i << "\n";
            std::cerr << "  expected " << primes << ", " << powers << "\n"
                      << "  got      " << primes_test << ", " << powers_test << "\n";
        }
    }
    return ok;
}

int main()
{
    std::cout << "test_example_1\n";
    test_example_1();
}


