#!/usr/bin/env bash

# path to executable
binary="build/bin/strategy_2/gen_nonrigid"

# maximum prime when constructing set of primes P
# ignored if primes are read in from file
max=100000000

#    min_size  - minimum size subsets of P to consider
#    max_size  - maximum size subsets of P to consider
#    L         - parameter L specified as a prime factorization
#               ex. "2^3 3^7" specifies L = 2^3 * 3^7
#    m1        - cofactor 1 (specifies the subgroups)
#    m2        - cofactor 2
#    m3        - cofactor 3

# usually maximum number of subsets to look at per partition (0 means no limit)
limit=1000

# the two non-rigid factors
p01="193 6337"

# minimum and maximum sizes of subets from each partition of P to consider
# (another way to limit the number of subsets besides $limit)
min_size=10
max_size=1000  # use 10000 or something large to include all

# factors of L
L_factors=" 2^5 3^5 5^4 7^2 11^1 13^1 17^1 19^1 29^1 31^1 37^1 41^1 47^1 53^1 59^1 61^1 67^1 71^1 73^1 79^1 83^1 89^1 97^1 101^1 "
# subgroups (cofactors m1 m2 and m3)
subgroups="7 11 13"

# path to read set of primes P from
# or '-' to generate the primes
primes="cache/3/primes"
# path to read set _a_ values from
# or '-' to generate them
a_values="cache/3/a_values"
 
$binary $max $limit $p01 $min_size $max_size "$L_factors"\
            $subgroups $primes $a_values

