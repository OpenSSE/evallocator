#!/usr/bin/env gnuplot --persist -c

if (ARGC < 1) {
    print "Usage: plot_load_n.gp <filename> or gnuplot plot_load_n.gp -c <filename>"
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

plot    "<(sed -e '/OneChoice/!d' ".filename.")" using 1:10 title "One choice mean maximum load" with linespoints,\
        "<(sed -e '/OneChoice/!d' ".filename.")" using 1:8 title "One choice max maximum load" with linespoints,\
        "<(sed -e '/OneChoice/!d' ".filename.")" using 1:9 title "One choice min maximum load" with linespoints, \
        "<(sed -e '/OneChoice/!d' ".filename.")" using 1:(3*$1/$2) title "One choice expected max load 3n/m" with lines,\
        "<(sed -e '/TwoChoice/!d' ".filename.")" using 1:10 title "Two choice mean maximum load" with linespoints,\
        "<(sed -e '/TwoChoice/!d' ".filename.")" using 1:8 title "Two choice max maximum load" with linespoints,\
        "<(sed -e '/TwoChoice/!d' ".filename.")" using 1:9 title "Two choice min maximum load" with linespoints, \
        "<(sed -e '/TwoChoice/!d' ".filename.")" using 1:(3*$1/$2) title "Two choice expected max load 4n/m" with lines
