Directory Structure
----------------
* `include`
   Contains all the header files for the project.

*  `plots`
   Contains collected data and generated plots from various runs.

*  `scripts`
   Collection of python and shell scripts for parsing data, scheduling jobs etc.

* `src`
   * `generate_cprimes.cpp`
   Program to generate Carmicheal number of any order using general Erdos construction.
   
   * `generate_cprimes_order_2.cpp`
   Program to generate Carmicheal number of order 2 using general Erdos construction.
   
   * `counting_factors.cpp`
   Module that implements methods that get the factorization of an integer using a sieve.
   
   * `gen_distributions.cpp`
   Prints out all possible exponents that could put on a list of provided numbers that still keeps the 
   magnitude of the product within a provided range.  
   
   * `timer.cpp`
   Module exports methods to allow for timing of sections of code.
   
   * `util.cpp`
   Collection of utility methods.
   
   * `construct_P.cpp`
   Implements two different algorithms for constructing a set of primes P for a given parameter L and a some limit.
   
   * `constuct_P_main.cpp`
   Program that parses command line arguments to get a value for the parameter L and a limit and uses the 
   specified algorithm from `construct_P.cpp` to construct a set of prime P.
   
   * `bench_construct_P.cpp`
   Program to benchmark the two algorithm impelementations in `construct_P.cpp` using the [google benchmarking library](https://github.com/google/benchmark).
 
 
 Building
----------------
All binaries can be built by running make with target `all`. By default all build files are placed in the generated `build` directory.

The following targets are currently defined in `Makefile`:

* `build/generate_cprimes`
   Target for building `generate_cprimes.cpp`.
   
* `build/generate_cprimes_order_2`
   Target for building `generate_cprimes_2.cpp`.

* `build/construct_P`
   Target for building `construct_P_main.cpp`.

* `build/bench_construct_P`
   Target for building `bench_construct_P.cpp`.
   
* `build/gen_distributions`
   Target for building `gen_distributions.cpp`.

The `build/` part is part of the target. For example, to build the benchmarking program:
```
$ make build/bench_construct_P
$ build/bench_construct_P
```
   
 

