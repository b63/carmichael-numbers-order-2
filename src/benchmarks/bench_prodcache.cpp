#include <vector>
#include <functional>
#include <cmath>

#include <NTL/ZZ.h>
#include <benchmark/benchmark.h>
#include <util.h>


void
subset_product_brute_force(const std::vector<long> &primes)
{
    size_t num_primes = primes.size();
    // go through all possible subset sizes starting from 2
    size_t count = 0;
    for (size_t t = 2; t <= num_primes; t++) 
    {
        size_t factors = t;
        std::vector<size_t> index_stack = {0};
#if LOG_LEVEL >= 1
        std::cout << "checking subsets of size " << t << " ...\n";
#endif

        /* cache of products */
        std::vector<NTL::ZZ> products; 
        products.reserve(t);

        // go through subsets of size factors suing a 'stack' (or odometer?)
        size_t index_size;
        while ((index_size = index_stack.size()) > 0) 
        {
#if LOG_LEVEL >= 2
            printVec<NTL::ZZ>(products);
            printVec<size_t>(index_stack);
            std::cout << "\n";
#endif
            const size_t prod_size  = products.size();
            size_t top = index_size - 1;
            size_t prod_top = prod_size - 1;

            if (index_stack[top] < num_primes) 
            {
                if (prod_size+1 < factors) 
                {
                    NTL::ZZ prod { 1 };
                    if (prod_size)
                        prod = products[prod_top];
                    products.push_back( prod * primes[index_stack[top]] );
                }

                if (index_size < factors) 
                {
                    index_stack.push_back(index_stack[top]+1);
                    continue;
                }
                /* IMPORTANT: must have at least TWO factors */
                /*            otherwise segfault  */
                NTL::ZZ prod = products[top - 1] * primes[index_stack[top]];
                benchmark::DoNotOptimize(prod);

# if LOG_LEVEL >= 1
                printVec<NTL::ZZ>(products);
                printVec<size_t>(index_stack);
                std::cout << " = " << prod << "\n";
#endif

                count++;
                index_stack[top]++;
            }
            else
            {
                index_stack.pop_back();
                if (top > 0) {
                    const size_t i = ++index_stack[top - 1];
                    products.pop_back();
                    if (i < num_primes) 
                    {
                        NTL::ZZ prod { 1 };
                        if (prod_size > 1)
                            prod = products[prod_top-1];
                        products.push_back(prod * NTL::ZZ{primes[i]});
                    }
                }

                continue;
            }
        }
    }

# if LOG_LEVEL >= 1
    std::cout << "total count: " << count << "\n";
#endif
}

void
subset_product_brute_force_no_cache(const std::vector<long> &primes)
{
    size_t num_primes = primes.size();
    // go through all possible subset sizes starting from 2
    size_t count = 0;
    for (size_t t = 2; t <= num_primes; t++) 
    {
        size_t factors = t;
        std::vector<size_t> index_stack = {0};
#if LOG_LEVEL >= 1
        std::cout << "checking subsets of size " << t << " ...\n";
#endif

        /* cache of products */
        std::vector<NTL::ZZ> products; 
        products.reserve(t);

        // go through subsets of size factors suing a 'stack' (or odometer?)
        size_t index_size;
        while ((index_size = index_stack.size()) > 0) 
        {
            size_t top = index_size - 1;

            if (index_stack[top] < num_primes) 
            {
                if (index_size < factors) 
                {
                    index_stack.push_back(index_stack[top]+1);
                    continue;
                }

                NTL::ZZ prod { 1 };
                for(size_t i = 0; i < index_size; i++)
                    prod *= primes[index_stack[i]];
                benchmark::DoNotOptimize(prod);

                count++;
                index_stack[top]++;
            }
            else
            {
                index_stack.pop_back();
                if (top > 0) {
                    ++index_stack[top - 1];
                }

                continue;
            }
        }
    }

# if LOG_LEVEL >= 1
    std::cout << "total count: " << count << "\n";
#endif
}

static void with_product_cache(benchmark::State &state)
{
    const size_t size = state.range(0);
    std::vector<long> set;
    set.reserve(size);
    for(size_t i = 0; i < size; i++)
        set.push_back((2L<<55) * (i+1));

    size_t subsets = pow(2, size)-size-1;
    for (auto _ : state)
    {
        subset_product_brute_force(set);
    }

    //printVec<NTL::ZZ>(products);
    std::cout << ".";
    state.counters["|P|"] = subsets;
}

static void without_product_cache(benchmark::State &state)
{
    const size_t size = state.range(0);
    std::vector<long> set;
    set.reserve(size);
    for(size_t i = 0; i < size; i++)
        set.push_back((2L<<55) * (i+1));

    size_t subsets = pow(2, size)-size-1;
    for (auto _ : state)
    {
        subset_product_brute_force_no_cache(set);
    }

    //printVec<NTL::ZZ>(products);
    std::cout << ".";
    state.counters["|P|"] = subsets;
}


BENCHMARK(with_product_cache)   ->DenseRange(4, 15, 1);
BENCHMARK(without_product_cache)->DenseRange(4, 15, 1);
