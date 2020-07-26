#include <iostream>
#include <vector>

#include <util.h>
#include <nonrigid.h>

int
main(int argc, char **argv)
{
    std::vector<std::vector<size_t> > cprimes;
    std::vector<long> primes {1, 2, 3, 4};
    subset_product_brute_force(cprimes, primes);
}
