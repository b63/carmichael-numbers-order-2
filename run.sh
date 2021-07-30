#!/usr/bin/env bash

# path to executable
binary="build/bin/strategy_2/gen_nonrigid"

# maximum prime when constructing set of primes P
# ignored if primes are read in from file
max=1000000000

# usually maximum number of subsets to look at per partition (0 means no limit)
limit=1000000

# the two non-rigid factors
p01="193 6337"

# minimum and maximum sizes of subets from each partition of P to consider
# (another way to limit the number of subsets besides $limit)
min_size=7
max_size=10  # use 10000 or something large to include all

# factors of L
L_factors="2^6 3^5 5^3 7^2 11 13 17 19 23 29 31 37 41 43 47 53"
# subgroups (cofactors m1 m2 and m3)
subgroups="18312281 14535931 266731860866625"

# path to read set of primes P from
# or '-' to generate the primes
primes="cache/7/primes"
# path to read set _a_ values from
# or '-' to generate them
a_values="cache/7/a_values"
 
$binary $max $limit $p01 $min_size $max_size "$L_factors"\
            $subgroups $primes $a_values

