#include <util.h>
#include <timer.h>
#include <NTL/ZZ.h>

int
test_example_1()
{
    long test_cases[5][2] {
        /*      n                 φ(n)         */
        {              1,                   1},
        {           2<<5,                  32},
        {             76,                  36},
        { 71234123412346,      35617061706172},
        {100000000000050,      26035409897280},
    };
    const size_t num_cases {sizeof(test_cases)/sizeof(test_cases[0])};

    int ok {1};
    for (size_t i {0}; i < num_cases; i++)
    {
        long *tcase {test_cases[i]};
        NTL::ZZ ret {eulers_toitent(NTL::ZZ{tcase[0]})};

        if (ret != tcase[1])
        {
            ok = 0;
            std::cerr << "FAIL case " << i << "\n";
            std::cerr << "  expected " << tcase[1] << ", got " << ret << "\n";
        }
    }
    return ok;
}


int
test_example_2()
{
    const char *cases[3] {
        "113133887486449488520580511944948209773520914929004851342263493591040000",
        "7651748237958879356666152743289756917852668142383579825413781282888888",
        "67521745127542765347625347523745237"
    };
    const NTL::ZZ answers[] {
        NTL::conv<NTL::ZZ>("13221612190310810376003999955754544436765178961223274264002560000000000"),
        NTL::conv<NTL::ZZ>("3659379368038969330444096868527029877342274640700163691127800118789632"),
        NTL::conv<NTL::ZZ>("63549877767089878325763558560277184"),
    };
    const size_t num_cases {sizeof(cases)/sizeof(cases[0])};

    int ok {1};
    init_timer();
    for(size_t i {0}; i < num_cases; ++i)
    {
        int id = start();
        NTL::ZZ ret { eulers_toitent( NTL::conv<NTL::ZZ>(cases[i]))};
        std::cout << end(id) << "\n";
        if (ret != answers[i])
        {
            ok = 0;
            std::cerr << "FAIL case " << i << "\n";
            std::cerr << "  expected " << answers[i] << ", got " << ret << "\n";
        }

    }
    return ok;
}


int main()
{
    std::cout << "test_example_1\n";
    test_example_1();
    test_example_2();
}


