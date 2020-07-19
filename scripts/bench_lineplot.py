import sys
import matplotlib.pyplot as plt
import matplotlib as m
import json

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

def cli_opt_log(ind, options):
    options['log_scale'] = True

def cli_opt_logbase(ind, options):
    options['log_scale_base'] = int(sys.argv[ind+1])

def cli_opt_output(ind, options):
    options['out_file'] = sys.argv[ind+1]

CLI_OPTIONS = {
        'l': [0, cli_opt_log],
        'o': [1, cli_opt_output],
        'b': [1, cli_opt_logbase],
    }
PLOT_VALUES = ['real_time']

def parse_cli_options():
    i = 1
    argc = len(sys.argv)
    options = {'out_file': 'benchmark.svg', 'log_scale': False, 'log_scale_base':2 }

    while i < argc:
        arg = sys.argv[i].strip()
        ind = arg.find('-')
        if ind < 0:
            break

        flag = arg[ind+1:]
        if flag in CLI_OPTIONS:
            num_args, func = CLI_OPTIONS[flag]
            if ind+num_args >= argc:
                not_enough_arguments(flag)

            func(i, options)
            i += num_args
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
        options['benchmark_file'] = sys.argv[i]

    return options

def main():
    options = parse_cli_options()
    bench_file = options['benchmark_file']
    log_scale = options['log_scale']
    image_output = options['out_file']

    # matplotlib variables
    time_units = None

    # collect data
    with open(bench_file, 'r') as f:
        benchjson =  json.loads(f.read())
    benchmarks = benchjson['benchmarks']

    # parse data
    benchmark_data = dict()
    for data in list(filter(lambda b: b['run_type'] == 'aggregate', benchmarks)):
        run_name, aggregate = data['run_name'], data['aggregate_name']
        parts = run_name.split('/')
        benchmark_name, run = parts[:2]

        try:
            run = int(run)
        except:
            print('Error: in "{}", "{}" is not a number'.format(run_name, run), file=sys.stderr)
            sys.exit(1)

        if benchmark_name not in benchmark_data:
            points = dict()
            benchmark_data[benchmark_name] = points
        else:
            points = benchmark_data[benchmark_name]

        if run not in points:
            #        mean, error
            point = [None, None]
            points[run] = point
        else:
            point = points[run]

        units = data['time_unit']
        if time_units is None:
            time_units = units
        elif time_units != units:
            print('Error: time units not the same for {}'.format(benchmark_name), file=sys.stderr)
            sys.exit(1)

        if aggregate == 'mean':
            point[0] = [data[v] for v in PLOT_VALUES]
        elif aggregate == 'stddev':
            point[1] = [data[v] for v in PLOT_VALUES]


    # plot the data
    fig, ax = plt.subplots()
    numcolors, cmap = 9, plt.get_cmap('Set1')

    colors = [cmap(i/numcolors) for i in range(9)]
    linestyles = ['solid', 'dashed']
    markers = ['.', 'o']
    num_plots = len(PLOT_VALUES)

    for i, v in enumerate(benchmark_data.items()):
        benchmark, points = v
        x      = []
        y_arr  = [[] for i in range(num_plots)]
        errors = [[] for i in range(num_plots)]

        for run, values in points.items():
            for j in range(num_plots):
                y_arr[j].append(values[0][j])
                errors[j].append(values[1][j])

            x.append(run)

        plot_lineargs = {'color': colors[i%numcolors], 'markersize':5}
        err_lineargs       = {'fmt':'none', 'color':'red', 'capsize':5, 'elinewidth':0.5, 'capthick':0.5}

        for j, plot in enumerate(PLOT_VALUES):
            if log_scale:
                ax.semilogx(x, y_arr[j], basex=options['log_scale_base'], 
                        label='{}({})'.format(benchmark, plot),
                        linestyle=linestyles[j%len(linestyles)],
                        marker=markers[j%len(markers)],
                        **plot_lineargs)
            else:
                ax.plot(x, y_arr[j], label='{}({})'.format(benchmark, plot),
                        marker=markers[j%len(markers)], 
                        linestyle=linestyles[j%len(linestyles)],
                        **plot_lineargs)
            ax.errorbar(x, y_arr[j], yerr=errors[j], **err_lineargs)

    ## increase range in y a litle bit
    #ymin, ymax = axes.get_ylim()
    #axes.set_ylim((ymin, ymax*1.5))

    lg = ax.legend(bbox_to_anchor=(1.05, 1))
    plt.title('Benchmarks')
    plt.ylabel('Time/max|P| ({})'.format(time_units))
    plt.xlabel('Limit')


    print('Saving figure to \'{}\'...'.format(image_output))
    plt.savefig(image_output, bbox_extra_artists=(lg,), bbox_inches='tight')

if __name__ == '__main__':
    main()







