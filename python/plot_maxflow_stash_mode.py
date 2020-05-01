#! python3

##
# Plot the experimental results of the buckets overflows
# The script takes a mandatory argument: the input file
# You can also specify the label to use when plotting by using the -l (or
# --label) option. The acceptable values are "m", "n", "n/m", "max_len" and
# "algorithm"

import argparse
import matplotlib.pyplot as plt
import json
import math


def plot_file(filename, label, normalize=False, logx=False, logy=False):
    with open(filename) as source:
        data = json.load(source)
        x = []
        y_max = []
        y_avg = []
        for experiment in data:  # iterate over the experiments
            n = float(experiment["parameters"]["exp_params"]["n"])
            m = float(experiment["parameters"]["exp_params"]["m"])
            p = int(experiment["parameters"]
                    ["exp_params"]["bucket_capacity"])

            iterations = int(experiment["parameters"]["iterations"])

            if label == "n/m":
                l = n / m
            elif label == "epsilon":
                l = (m*p/n) - 2
            else:
                l = float(experiment["parameters"]["exp_params"][label])

            norm_factor = 1

            if (normalize):
                norm_factor = n

            plt.plot([i/(norm_factor*iterations) for i in experiment["stash_modes"]
                      [p:: p]], label="n=%d" % (n))
            # x.append(l)

            # stash_max = float(experiment["stash_size"]["max"])
            # stash_avg = float(experiment["stash_size"]["mean"])

            # if normalize:
            # stash_avg /= n
            # stash_max /= n

            # y_max.append(stash_max)
            # y_avg.append(stash_avg)

        # plt.plot(x, y_max, label="Max stash size")
        # plt.plot(x, y_avg, label="Average stash size")
        plt.legend(loc='upper right')

        if logx:
            plt.semilogx()
        if logy:
            plt.semilogy()

        plt.show()


parser = argparse.ArgumentParser(description='Plot allocation overflows.')
parser.add_argument('filename', metavar='path',
                    help='Path to a JSON file')
parser.add_argument('--label', '-l', default='n',
                    help='Define the used label')
parser.add_argument('--normalize', '-n', action='store_true')
parser.add_argument('--logx', action='store_true')
parser.add_argument('--logy', action='store_true')


args = parser.parse_args()
# print(args)

plot_file(args.filename, args.label, normalize=args.normalize,
          logx=args.logx, logy=args.logy)
