import matplotlib.pyplot as plt
import matplotlib as m

def interp(a, b, f):
    mid = [ a[i]+f*(b[i]-a[i]) for i in range(len(a))]
    return mid


SUMMARY_FILE='data/inter/summary'
MAX=200

data = []
primes = []
with open(SUMMARY_FILE, 'r') as f:
    i = 0
    while True:
        line = f.readline().strip()
        if len(line) < 1:
            break
        elif (i == 0):
            for s in line.split(' '):
                if len(s) > 0:
                    primes.append(int(s))
            i += 1
            continue

        values = line.split(',')
        if len(values) != 3:
            print('ignoring line, {}'.format(line))
            continue
        dist = [ int(n) for n in values[1].strip().split(' ')]
        density = float(values[2])
        data.append((density, dist))
        i+=1

data.sort(key=lambda v: v[0], reverse=True)

# plot data
cm = m.cm.get_cmap('RdPu')
nm = plt.Normalize(data[min(MAX,len(data)-1)][0], data[0][0], clip=True)
sm = m.cm.ScalarMappable(cmap=cm, norm=nm)

x = [i for i in range(len(data[0][1]))]
y2 = [0 for i in range(len(data[0][1]))]

alpha = 1
factor = 0.5
for i in range(len(data)):
    if alpha < 0.1 or i >= MAX:
        break

    density, points = data[i]
    clr = list(cm(nm(density)))
    clr[-1] = alpha

    print(points)
    plt.plot(x, points, alpha=alpha,color=clr, zorder=alpha)
    alpha -= factor/(i+1)
    factor = 0.2

cbar = plt.colorbar(sm)
cbar.set_label('Density')

plt.ylabel('Power');
plt.xlabel('Index of Prime number');
plt.xticks([i for i in range(len(primes))], [str(p) for p in primes])
if MAX < len(data):
    plt.title('{} Distributions with Highest density'.format(MAX))
else:
    plt.title('Distributions ({})'.format(len(data)))
plt.savefig('test.png')
