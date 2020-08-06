#include <util.h>

int
test_example_1()
{
    for(size_t j {0}; j < 10; ++j)
    {
        std::unique_ptr<std::vector<size_t>> ind {random_subset(10)};
        std::cout << j << ": " << *ind << "\n";
    }

    return 1;
}

int main()
{
    std::cout << "test_example_1\n";
    test_example_1();
}


