#!/usr/bin/env gnuplot --persist -c

if (ARGC < 1) {
    print "Usage: plot_load_n_m.gp <filename> or gnuplot plot_load_n_m.gp -c <filename>"
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

# set logscale x 2


# Uncomment the following lines to run a linear regression
# f(x) = a*x + b
# fit f(x) "<(sed -e '/,TwoChoice/!d' ".filename.")" using ($1/$2):10 via a, b


plot    "<(sed -e '/,OneChoice/!d' ".filename.")" using ($1/$2):12 title "One choice mean maximum load" with linespoints,\
        "<(sed -e '/,OneChoice/!d' ".filename.")" using ($1/$2):10 title "One choice max maximum load" with linespoints,\
        "<(sed -e '/,OneChoice/!d' ".filename.")" using ($1/$2):11 title "One choice min maximum load" with linespoints, \
        "<(sed -e '/,OneChoice/!d' ".filename.")" using ($1/$2):(3*$1/$2) title "One choice expected max load 3n/m" with lines,\
        "<(sed -e '/,BlockedOneChoice/!d' ".filename.")" using ($1/$2):12 title "Blocked One choice mean maximum load" with linespoints,\
        "<(sed -e '/,BlockedOneChoice/!d' ".filename.")" using ($1/$2):10 title "Blocked One choice max maximum load" with linespoints,\
        "<(sed -e '/,BlockedOneChoice/!d' ".filename.")" using ($1/$2):11 title "Blocked One choice min maximum load" with linespoints, \
        "<(sed -e '/,TwoChoice/!d' ".filename.")" using ($1/$2):12 title "Two choice mean maximum load" with linespoints,\
        "<(sed -e '/,TwoChoice/!d' ".filename.")" using ($1/$2):10 title "Two choice max maximum load" with linespoints,\
        "<(sed -e '/,TwoChoice/!d' ".filename.")" using ($1/$2):11 title "Two choice min maximum load" with linespoints




# additional graphs :
# "<(sed -e '/,OneChoice/!d' ".filename.")" using ($1/$2):10:8:9 notitle with errorbars,\ # error bars based on min and max
# "<(sed -e '/,OneChoice/!d' ".filename.")" using ($1/$2):10:1($1/$2):9 notitle with errorbars,\ # error bars based on std dev
# "<(sed -e '/,TwoChoice/!d' ".filename.")" using ($1/$2):(f($1/$2)) title "fit" with lines,\ # print the fitted function