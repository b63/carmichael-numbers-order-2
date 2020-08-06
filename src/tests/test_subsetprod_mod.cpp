
#include <util.h>
#include <strategy_2/subset_product.h>


template <unsigned long int N>
int
check_set_base(const std::vector<long> &set, const std::array<NTL::ZZ, N> &bases, int verbose)
{
    constexpr size_t num_bases {N};

    int retv {1};
    subsetprod_mod<num_bases>(set, bases, 
            [&](std::array<NTL::ZZ, num_bases>& prod, const std::vector<size_t>& indices, size_t insert_index) -> bool {
                NTL::ZZ p {1};
                for(auto it {indices.cbegin()};it!=indices.cend(); ++it)
                {
                    p *= set[*it];
                    if (verbose)
                    {
                        std::cout << set[*it];
                        if (it+1 != indices.cend())
                            std::cout << " * ";
                    }
                }

                if (verbose)
                    std::cout << " = " << p << " " << prod
                            << " mod " << bases << " ";
                NTL::ZZ r;
                int ok {1};
                for (size_t i {0}; i < num_bases; ++i)
                {
                    NTL::rem(r, p, bases[i]);
                    if(verbose)
                        std::cout << ".";
                    if (r != prod[i])
                    {
                        retv = ok = 0;
                        std::cout << "FAIL\n  ";
                        std::cout << "indices: " << indices << ", mod: " << bases[i] 
                            << ", expected " << r << " got " << prod[i] << "\n";
                        break;
                    }
                }
                if(verbose && ok)
                    std::cout << " PASS\n";
            return false;
            }, 1, 0);
    return 1;
}

int
test_example_1()
{
    std::vector<long> set {2, 3, 4, 5, 7};
    std::array<NTL::ZZ, 3> bases {
        NTL::ZZ{2},
        NTL::ZZ{4},
        NTL::ZZ{9}
    };

    return check_set_base<3>(set, bases, 0);
}


int
test_example_2()
{
    std::vector<long> set {2, 3, 4, 5, 7, 12, 45, 223144, 1238941876871234L, 187129389823};
    std::array<NTL::ZZ, 5> bases {
        NTL::ZZ{293847192384192834},
        NTL::conv<NTL::ZZ>("9812934918891234981723948172839412734"),
        NTL::ZZ{21341241234},
        NTL::ZZ{241234},
        NTL::ZZ{24988912391234},
    };

    return check_set_base<5>(set, bases, 0);
}


int
test_example_3()
{
    std::vector<long> set {187129389823};
    std::array<NTL::ZZ, 2> bases {
        NTL::ZZ{2},
        NTL::ZZ{3},
    };

    return check_set_base<2>(set, bases, 0);
}

int main()
{
    std::cout << "test_example_1\n";
    test_example_1();
    std::cout << "test_example_2\n";
    test_example_2();
    std::cout << "test_example_3\n";
    test_example_2();
}


