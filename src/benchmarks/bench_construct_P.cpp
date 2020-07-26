#include <vector>
#include <iomanip>

#include <NTL/ZZ.h>
#include <benchmark/benchmark.h>

#include <construct_P.h>
#include <util.h>
#include <timer.h>

static std::vector<long> L_P_primes        = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
static std::vector<long> L_P_primes_powers = {15, 6, 3, 1,  1,  1,  1,  1,  1,  1,  1,  1 };

static void Method_2(benchmark::State &state)
{
    size_t lim = state.range(0);
    init_timer();

    // multiply the prime factors with appropriate powers to get L'
    NTL::ZZ L_P { 1 };
    multiply_factors(L_P, L_P_primes, L_P_primes_powers);
    size_t k_max_size = 0;

    for (auto _ : state) 
    {
        std::map<const long, std::vector<NTL::ZZ> > factor_map{};
        lim = lim*lim;
        int id { start(-1) };
        construct_primes_2(factor_map, L_P_primes, L_P_primes_powers, lim);
        time_metric times { end(id) };
        benchmark::DoNotOptimize(factor_map);

        // find max size
        size_t max_size = 0;
        for (auto it=factor_map.cbegin(); it != factor_map.cend(); it++)
        {
            size_t size = it->second.size();
            if (size > max_size)
                max_size = size;
        }

        k_max_size = max_size;
        std::cout << std::setprecision(10) << times.wall_time/max_size << "\n";
        state.SetIterationTime(times.wall_time/max_size);
    }
    state.counters["|P|"] = k_max_size;
}

static void Method_1(benchmark::State &state)
{
    size_t lim = state.range(0) ;
    init_timer();

    // multiply the prime factors with appropriate powers to get L'
    NTL::ZZ L_P { 1 };
    multiply_factors(L_P, L_P_primes, L_P_primes_powers);
    size_t k_max_size = 0;

    for (auto _ : state) 
    {
        std::map<const long, std::vector<long> > factor_map{};

        int id { start() };
        construct_primes(factor_map, L_P_primes, L_P_primes_powers, lim);
        time_metric times { end(id) };

        benchmark::DoNotOptimize(factor_map);

        // find max size
        size_t max_size = 0;
        for (auto it=factor_map.cbegin(); it != factor_map.cend(); it++)
        {
            size_t size = it->second.size();
            if (size > max_size)
                max_size = size;
        }

        k_max_size = max_size;
        std::cout << std::setprecision(10) << times.wall_time/max_size << "\n";
        state.SetIterationTime(times.wall_time/max_size);
    }
    state.counters["|P|"] = k_max_size;
}

BENCHMARK(Method_1)->Unit(benchmark::kMillisecond)->Range(1024, 1024<<15)->UseManualTime();
BENCHMARK(Method_2)->Unit(benchmark::kMillisecond)->Range(1024, 1024<<15)->UseManualTime();
