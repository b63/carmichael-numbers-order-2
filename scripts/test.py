import sys
from functools import reduce

f = sys.stdin

line = f.readline()
factors = [int(x) for x in line.split(' ')]
n = reduce(lambda x,y: x*y, factors)
print('n = {:d}'.format(n))
for p in factors:
    print("p = {:d}, n mod p^2-1 = {}".format(p, (n)%(p**2-1)))


