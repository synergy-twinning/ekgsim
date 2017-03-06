#!/bin/bash

# get terminal size
temp="$(mktemp)"
resize > $temp
source $temp
cat $temp
rm $temp

./hyper.sh ${1-"individuals.txt"} | gnuplot -e "set term dumb $COLUMNS $LINES; set xlabel 'generations'; set ylabel 'hypervolume'; set key off; plot '-' using 1:2 w lines"
