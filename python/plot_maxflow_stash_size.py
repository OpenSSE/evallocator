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


def plot_file(filename, label):
    with open(filename) as source:
        data = json.load(source)
        x = []
        y_max = []
        y_avg = []
        for experiment in data:  # iterate over the experiments
            n = float(experiment["parameters"]["exp_params"]["n"])
            m = float(experiment["parameters"]["exp_params"]["m"])
            p = float(experiment["parameters"]
                      ["exp_params"]["bucket_capacity"])

            if label == "n/m":
                l = n / m
            elif label == "epsilon":
                l = (m*p/n) - 2
            else:
                l = float(experiment["parameters"]["exp_params"][label])

            x.append(l)

            y_max.append(float(experiment["stash_size"]["max"]))
            y_avg.append(float(experiment["stash_size"]["mean"]))

        plt.plot(x, y_max, label="Max stash size")
        plt.plot(x, y_avg, label="Average stash size")
        plt.legend(loc='upper left')
        plt.show()


parser = argparse.ArgumentParser(description='Plot allocation overflows.')
parser.add_argument('filename', metavar='path',
                    help='Path to a JSON file')
parser.add_argument('--label', '-l', default='n',
                    help='Define the used label')

args = parser.parse_args()
# print(args)

# filename = "../experiments/var_n_m/large_m_res.json"
# filename = "../experiments/var_max_len/one_choice_res.json"

plot_file(args.filename, args.label)
