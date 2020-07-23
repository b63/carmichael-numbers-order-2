#include <iostream>
#include <vector>

#include <util.h>
#include <nonrigid.h>

int
main(int argc, char **argv)
{
    Factorization f1 { {2, 3, 5, 7}, {2, 2, 2, 2} };
    Factorization f2 { {2, 3, 5, 7}, {2, 2, 2, 2} };

    std::cout << f1 << "\n" << f2 << "\n";
    std::cout << include_as_factor(f1, f2) << "\n";
}
