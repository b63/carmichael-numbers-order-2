import sys
from functools import reduce

f = sys.stdin

line = f.readline()
factors = [int(x) for x in line.split(' ')]
n = reduce(lambda x,y: x*y, factors)
for p in factors:
    print("p = {:d}, n mod p^2-1 = {}".format(p, (n-1)%(p**2-1)))


