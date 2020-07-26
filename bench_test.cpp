#include <vector>
#include <benchmark/benchmark.h>


static void Without_Reserve(benchmark::State &state)
{

    for (auto _ : state)
    {
        std::vector<long> data;
        for (size_t i = 0; i < 100; i++)
            data.push_back(100);
        benchmark::DoNotOptimize(data);
    }
}

static void With_Reserve(benchmark::State &state)
{

    for (auto _ : state)
    {
        std::vector<long> data;
        data.reserve(100);
        for (size_t i = 0; i < 100; i++)
            data.push_back(100);
        benchmark::DoNotOptimize(data);
    }
}

BENCHMARK(With_Reserve);
BENCHMARK(Without_Reserve);
