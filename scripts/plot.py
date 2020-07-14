import sys
import matplotlib.pyplot as plt
import matplotlib as m

def interp(a, b, f):
    mid = [ a[i]+f*(b[i]-a[i]) for i in range(len(a))]
    return mid

if len(sys.argv) < 2:
    print("Please provide path to data summary file (data/interpolation/summary)")
    sys.exit(1)

SUMMARY_FILE=sys.argv[1]
OUT_FILE='plot.svg'
SORTED_FILE='sorted_summary'
PRIME_LABEL=False
MAX=50

data = []
primes = []
options=dict()
with open(SUMMARY_FILE, 'r') as f:
    i = 0
    while True:
        line = f.readline()
        if len(line) < 1:
            break
        line = line.strip()
        if (i == 0):
            for s in line.split(' '):
                if len(s) > 0:
                    primes.append(int(s))
            i += 1
            continue
        elif (i==1):
            for s in line.split(' '):
                s  = s.strip()
                if len(s) < 2: continue
                flag = s[:2]
                if flag == '-m':
                    options['Method'] = s[2:]
                elif flag == '-l':
                    options['Limit'] = s[2:]

        values = line.split(',')
        if len(values) != 3:
            print('ignoring line, {}'.format(line))
            continue
        dist = [ int(n) for n in values[1].strip().split(' ')]
        density = float(values[2])
        data.append((density, dist))
        i+=1

data.sort(key=lambda v: v[0], reverse=True)

# write sorted data to file
with open(SORTED_FILE, 'w') as f:
    print('Saving sorted data to \'{}\'...'.format(SORTED_FILE))
    for i, v in enumerate(data):
        dist = ' '.join([str(_) for _ in v[1]])
        f.write('{:d}, {:.10f}, {}\n'.format(i, v[0], dist))

# plot data
cm = m.cm.get_cmap('RdPu')
nm = plt.Normalize(data[min(MAX,len(data)-1)][0], data[0][0], clip=True)
sm = m.cm.ScalarMappable(cmap=cm, norm=nm)

x = [i for i in range(len(data[0][1]))]
y2 = [0 for i in range(len(data[0][1]))]

alpha, minalpha = 1, 0.1
factor = 0.5
for i in range(len(data)):
    if alpha < 0.01 or i >= MAX:
        break

    density, points = data[i]
    clr = list(cm(nm(density)))
    clr[-1] = alpha

    plt.plot(x, points, alpha=alpha,color=clr, zorder=alpha)
    alpha= max(alpha - factor/(i+1), minalpha)
    factor = 0.2

cbar = plt.colorbar(sm)
cbar.set_label('Density')

plt.grid(True, linestyle='dotted', alpha=0.5)
plt.ylabel('Power');
plt.axis(xmin=0, ymin=0, xmax=len(primes)-1)
ax, fig = plt.gca(), plt.gcf()

if PRIME_LABEL:
    plt.xlabel('Prime number');
    plt.xticks([i for i in range(len(primes))], [str(p) for p in primes])
else:
    plt.xlabel('Index');
    ax.set_xticks(ticks=[2*i for i in range(len(primes)//2)], minor=False)
    ax.set_xticks(ticks=[i for i in range(len(primes))], minor=True)
    plt.subplots_adjust(bottom=0.2)

    txt = 'Primes: [' + ', '.join([str(p) for p in primes]) + '] '
    txt += '\n' + ', '.join([k+'='+str(v) for k, v in options.items()])
    textobj = plt.figtext(0.1, 0.03, txt, transform=plt.gcf().transFigure)

    fig = plt.gcf()
    fig.canvas.draw()
    txt_bbox = textobj.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    if fig.get_figwidth() < txt_bbox.xmax:
        fig.set_figwidth(txt_bbox.xmax + 0.5)


if MAX < len(data):
    plt.title('{} Distributions with Highest density'.format(MAX))
else:
    plt.title('Distributions ({})'.format(len(data)))
print('Saving figure to \'{}\'...'.format(OUT_FILE))
plt.savefig(OUT_FILE)
