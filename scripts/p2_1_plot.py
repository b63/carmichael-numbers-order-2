import sys
import matplotlib.pyplot as plt
import matplotlib as m
import json

def increase_axes_y(axes):
    # increase range in y a litle bit
    ymin, ymax = axes.get_ylim()
    axes.set_ylim((ymin, ymax*1.5))

INPUT_FILE=None
OUT_FILE='p2_1.svg'


# parse CLI arguments
argc = len(sys.argv)
if argc > 2:
    print('Error: too many arguments', file=sys.stderr)
    sys.exit(1)
elif argc < 2:
    print('Error: no input file specified', file=sys.stderr)
    sys.exit(1)
else:
    INPUT_FILE = sys.argv[1]


# matplotlib variables
bar_width = 1.5

# collect data
pdata = dict()
primes = []
num_points = 0
with open(INPUT_FILE, 'r') as f:
    readprimes = False
    while (line:=f.readline()):
        if not readprimes and line.startswith('primes'):
            readprimes = True
            for s in line[line.find(':')+1:].split(' '):
                s = s.strip()
                if s.isnumeric(): 
                    primes.append(int(s))
            print('read primes: {}'.format(primes))
            continue

        ind = line.find('=')
        if ind == -1:
            print('ignoring line: {} ...'.format(line[:10]))
            continue
        line = line[ind+1:]
        factors = line.split(' ')
        for factor in factors:
            arr = factor.split('^')
            if len(arr) != 2:
                continue
            prime, power = [int(i) for i in arr]
            if prime in pdata:
                data = pdata[prime]
                data[0] += power
                data[1] = max(power, data[1])
                data[2] = min(power, data[2])
            else:
                data = [power, power, power]
                pdata[prime] = data
        num_points += 1

print('parsed {} lines ...'.format(num_points))

# put data into correct format for matplotlib
x = [i for i in range(1, len(primes)+1)]
y_avg, y_max = [], []
for p in primes:
    _avg, _max, _min = pdata.get(p, [0,0,0])
    y_avg.append(_avg/num_points)
    y_max.append(_max)

#print(list(zip(x, y_avg, y_max)))

title = 'Factorization of $p^2-1$ for first {} primes'.format(num_points)

# increase figure size
fig = plt.gcf()
w, h = fig.get_size_inches()
fig.set_size_inches(w*max(1, num_points/120), h)

# plot average data
bars = plt.bar(x, y_avg, align='center', width=0.5)
increase_axes_y(plt.gca())
plt.grid(True, linestyle='dotted', alpha=0.5)
plt.title(title)
plt.ylabel('Average power')
plt.xlabel('Rank of Prime')

outfile = OUT_FILE.replace('.', '_avg.')
print('Saving figure to \'{}\'...'.format(outfile))
plt.savefig(outfile)

# plot max data
plt.cla()
bars = plt.bar(x, y_max, align='center', width=0.5)
increase_axes_y(plt.gca())
plt.grid(True, linestyle='dotted', alpha=0.5)
plt.title(title)
plt.ylabel('Max power')
plt.xlabel('Rank of Prime')

outfile = OUT_FILE.replace('.', '_max.')
print('Saving figure to \'{}\'...'.format(outfile))
plt.savefig(outfile)


