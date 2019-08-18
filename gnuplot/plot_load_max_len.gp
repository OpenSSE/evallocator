#!/usr/bin/env gnuplot --persist -c

if (ARGC < 1) {
    print "Usage: plot_load_max_lean.gp <filename> or gnuplot plot_load_max_lean.gp -c <filename>"
    exit
}

filename=ARG1

file_exists(file) = system("[ -f '".file."' ] && echo '1' || echo '0'") + 0

if (!file_exists(filename)) {
    print "Invalid path. File \"".filename."\" does not exist"
    exit
}

reset
set datafile separator ","

set logscale x 2

plot    filename using 3:10 title "Mean maximum load" with linespoints,\
        filename using 3:8 title "Max maximum load" with linespoints,\
        filename using 3:9 title "Min maximum load" with linespoints