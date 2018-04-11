#!/bin/bash

# get terminal size
temp="$(mktemp)"
resize > $temp
source $temp
cat $temp
rm $temp

# parse last front and draw it
./front.sh ${1:-"individuals.txt"} ${2:-"-1"} | gnuplot -e "set term dumb $COLUMNS $LINES; set xrange [0:*]; set yrange [0:*]; set xlabel 'objective 1'; set ylabel 'objective 2'; set key off; plot '-' using 1:2 w points"
