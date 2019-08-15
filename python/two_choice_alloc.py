import random
import numpy

from progressbar import ProgressBar, Counter, Timer, Bar, Percentage, AdaptiveETA


def power_of_two(x):
    # returns the smallest power of two larger than x
    return 1 if x == 0 else 2**(x - 1).bit_length()


def alloc(N, m, max_size):
    assert(m > max_size)
    # m is supposed to be a power of 2
    # create an array of m buckets
    # because we are only interested in the size of the buckets
    # there is no need to do something different from an array of integers
    buckets = numpy.zeros(m, numpy.int64)
    remaining_elements = N

    widgets = ['Processed: ', Counter(), ' elements (', Percentage(), ')',
               Bar(), Timer(), AdaptiveETA()]
    pbar = ProgressBar(widgets=widgets, maxval=N)
    pbar.start()

    while remaining_elements != 0:
        # generate a random list length
        l = random.randint(1, min(max_size, remaining_elements))

        # l = int(random.normalvariate(max_size/2, 3))
        # while (1 > l or l > min(max_size, remaining_elements)):
        # l = int(random.normalvariate(max_size / 2, 3))

        # print(l)
        # generate the two locations
        # remember that the list length are supposed to be power of 2's
        # we will not pad lists here, but the 'meta buckets' must still have a
        # size that is a power of 2
        n_i = power_of_two(l)
        meta_buckets_counts = m / n_i

        B_1 = random.randint(0, meta_buckets_counts-1)
        B_2 = random.randint(0, meta_buckets_counts-1)

        count_B_1 = sum(buckets[n_i*B_1: n_i*(B_1+1)])
        count_B_2 = sum(buckets[n_i*B_2: n_i*(B_2+1)])

        # choose the lesser full
        chosen_B = B_1
        if (count_B_1 > count_B_2):
            chosen_B = B_2

        for j in range(n_i * chosen_B, n_i * chosen_B + n_i):
            buckets[j] += 1

        remaining_elements -= l
        pbar.update(N-remaining_elements)

    return buckets


def print_experiment(N, m, max_len):
    allocation = alloc(N, m, max_len)
    tot = sum(allocation)
    max_load = max(allocation)

    print(tot)
    print(tot/N)
    print(max_load)
    print(allocation)


def iterated_experiment(N, m, max_len, iterations):
    max_load = 0

    print("Starting experiments...\n")

    for i in range(iterations-1):
        allocation = alloc(N, m, max_len)
        max_load = max(max_load, max(allocation))
        print(i + "\r")
    return max_load


if __name__ == "__main__":
    # N = 2 ** 10
    # m = 2 ** 5
    # max_len = 2**4

    N = 2 ** 30
    m = 2 ** 20
    max_len = 2**17

    print_experiment(N, m, max_len)
    # max = iterated_experiment(N, m, max_len, 1000)

    # print(max)
