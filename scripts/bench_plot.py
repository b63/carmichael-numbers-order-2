import sys
import matplotlib.pyplot as plt
import matplotlib as m
import json


BENCH_FILE=None
OUT_FILE=None
RELATIVE=False


def not_enough_arguments(flag):
    print('Error: no value provided for flag \'{}\''.format(flag), file=sys.stderr)
    sys.exit(1)

def invalid_flag(arg):
    print('Error: invalid option \'{}\''.format(arg), file=sys.stderr)
    sys.exit(1)

def label_bars(axes, bars, offset=1, suffix=''):
    """
    Label bars in barplot
    """
    for bar in bars:
        y = bar.get_height()
        x = bar.get_x() + bar.get_width()/2
        axes.text(x, y + offset, '{:.1f}{}'.format(y, suffix), ha='center', color='black', fontsize=7)


# parse CLI arguments
i = 1
argc = len(sys.argv)
while i < argc:
    arg = sys.argv[i].strip()
    ind = arg.find('-')
    flag = None
    if ind < 0:
        break
    else:
        flag = arg[ind+1:]

    if flag == 'r':
        RELATIVE = True
    elif flag == 'o':
        if i+1<argc:
            OUT_FILE = sys.argv[i+1]
            i += 1
        else:
            not_enough_arguments(arg)
    else:
        invalid_flag(arg)

    i += 1

if i+1 < argc:
    print('Error: too many arguments', file=sys.stderr)
    sys.exit(1)
elif i+1 > argc:
    print('Error: no benchmark json file specified', file=sys.stderr)
    sys.exit(1)
else:
    BENCH_FILE = sys.argv[i]

if OUT_FILE is None:
    OUT_FILE = 'benchmark.svg'



# matplotlib variables
wall_color, cpu_color, median_color, stddev_color = '#65C280', '#FBD650', '#002953', '#f00'
bar_width = 1.5
time_units = None

# collect data
with open(BENCH_FILE, 'r') as f:
    benchjson =  json.loads(f.read())
benchmarks = benchjson['benchmarks']

# parse data
unique_points = 0
points = dict()
reference, reference_bench=None, None
for data in list(filter(lambda b: b['run_type'] == 'aggregate', benchmarks)):
    benchmark_name, aggregate = data['run_name'], data['aggregate_name']
    if benchmark_name not in points:
        point = dict()
        points[benchmark_name] = point
    else:
        point = points[benchmark_name]

    units = data['time_unit']
    if time_units is None:
        time_units = units
    elif time_units != units:
        print('Error: time units not the same for {}'.format(benchmark_name), file=sys.stderr)
        sys.exit(1)


    if aggregate == 'mean':
        xstart = unique_points * 5
        point['x'] = (xstart, xstart + 2)
        point['height'] = (data['real_time'], data['cpu_time'])
        point['color'] = [wall_color, cpu_color]
        point['tick'] = (xstart+1, benchmark_name)

        # update reference
        if reference is None:
            reference = (data['real_time'], data['cpu_time'])
            reference_bench = benchmark_name

        unique_points += 1
    elif aggregate == 'stddev':
        point['yerr'] = (data['real_time'], data['cpu_time'])
    elif aggregate == 'median':
        point['median'] = (data['real_time'], data['cpu_time'])

# put data into correct format for matplotlib
x, height, yerr, = [], [], []
colors, x_ticks, labels = [], [], []
median_y, median_xmin, median_xmax = [], [], []
for name, point in points.items():
    x.extend(point['x'])

    heights = point['height']
    yerrs   = point['yerr']
    medians = point['median']

    if RELATIVE:
        heights = [heights[i]/reference[i] for i in range(len(reference))]
        medians = [medians[i]/reference[i] for i in range(len(reference))]
        yerrs   = [  yerrs[i]/reference[i] for i in range(len(reference))]

    height.extend(heights)
    yerr.extend(yerrs)
    median_y.extend(medians)

    median_xmin.extend([v-bar_width/2 - 0.2 for v in point['x']])
    median_xmax.extend([v+bar_width/2 + 0.2 for v in point['x']])

    colors.extend(point['color'])

    xtick, label = point['tick']
    x_ticks.append(xtick)
    labels.append(label)


# plot data
wall_patch  = m.patches.Patch(color=wall_color, label='Real Time')
cpu_patch   = m.patches.Patch(color=cpu_color, label='CPU Time')
median_patch = m.lines.Line2D([], [], color=median_color, linestyle='dotted', label='Median')
stddev_patch = m.lines.Line2D([], [], color=stddev_color, label='Std. dev')

bars = plt.bar(x, height, align='center', width=1.5, color=colors)
plt.hlines(median_y, median_xmin, median_xmax, colors=median_color, linestyles='dotted')
plt.xticks(x_ticks, labels)

axes = plt.gca()
axes.errorbar(x, height, yerr=yerr, fmt='none', ecolor=stddev_color, capsize=2)

# increase range in y a litle bit
ymin, ymax = axes.get_ylim()
axes.set_ylim((ymin, ymax*1.5))

# plot bar height labels
if RELATIVE:
    label_bars(axes, bars, offset=0.5, suffix='x')

plt.legend(handles=[wall_patch, cpu_patch, median_patch, stddev_patch])
plt.title('Benchmarks')
if RELATIVE:
    plt.ylabel('Time (Relative to {})'.format(reference_bench))
else:
    plt.ylabel('Time ({})'.format(time_units))



print('Saving figure to \'{}\'...'.format(OUT_FILE))
plt.savefig(OUT_FILE)







