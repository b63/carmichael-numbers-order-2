#include <strategy_2/subset_product.h>

int
test_example_1()
{
    size_t test_cases[5][4] {
        /* set_size, min_size, max_size, expected result */
        {10,  1,  1,                 10},
        {10,  1, 10,               1023},
        {10, 10, 10,                  1},
        {20, 14, 15,              54264},
        {50,  1, 50,   1125899906842623},
    };

    int ok {1};
    for (size_t i {0}, len=sizeof(test_cases)/sizeof(test_cases[0]); i < len; i++)
    {
        size_t *tcase {test_cases[i]};
        size_t ret {calc_max_subsets(tcase[0], tcase[1], tcase[2])};

        if (ret != tcase[3])
        {
            ok = 0;
            std::cerr << "FAIL case " << i << "\n";
            std::cerr << "  expected " << tcase[3] << ", got " << ret << "\n";
        }
    }
    return ok;
}

int main()
{
    std::cout << "test_example_1\n";
    test_example_1();
}


