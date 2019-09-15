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
    cmap = plt.get_cmap('jet')
    N = 20

    fig, (ax1) = plt.subplots(nrows=1, ncols=1)
    # fig = plt.figure()
    # ax1 = fig.add_subplot(211)
    # ax2 = fig.add_subplot(212)

    with open(filename) as souce:
        data = json.load(souce)
        exp_num = 10
        for experiment in data:  # iterate over the experiments
            stat_min = []
            stat_max = []
            stat_mean = []
            stat_std_dev = []
            count = 0
            # iterate over the different modes
            for mode in experiment["load_modes"]:
                # if (overflow_stat[1] == 0):
                    # break
                count = count+1
                stat_min.append(mode[0])
                stat_max.append(mode[1])
                stat_mean.append(mode[2])
                stat_std_dev.append(math.sqrt(mode[3]))

            color = cmap(float(exp_num)/N)
            # plt.plot(stat_min, dashes=[6, 2], color=color)
            l = ""
            if label == "n/m":
                l = str(float(experiment["parameters"]["n"]) /
                        float(experiment["parameters"]["m"]))
            else:
                l = experiment["parameters"][label]
            ax1.plot(range(0, count), stat_mean,
                     label=l, color=color, marker="x")
            exp_num = exp_num + 1

    # ax1.semilogy()
    # ax2.semilogy()
    # plt.xlabel('alpha')

    # plt.subplot(211)
    # plt.semilogy()
    # # plt.xlabel('alpha')
    # plt.ylabel('# overflowing entries')
    # # plt.legend()
    # plt.subplot(212)
    # plt.semilogy()
    # plt.xlabel('alpha')
    # ax1.set_ylabel('Max # overflows')
    # ax2.set_ylabel('Average # overflows')
    # ax3.set_ylabel('Max # overflows')
    # ax4.set_ylabel('Average # overflows')
    plt.legend()

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
