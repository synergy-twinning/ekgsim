#!/bin/bash

# takes individuals file as an input and outputs the last front (pace separated riteria) 
# awk removes the number of generation (first field in the return of AMS-DEMO) and comments (lines starting with #)
# use: extract the last front for plotting with matlab

./AMS-DEMO/DEMO -individuals:${1-"individuals.txt"} -analysis:front -gen:${2-"-1"} -out:criteria | awk '{s=""; for (i=2; i<=NF; i++) {s=s $i; if (i!=NF) s=s " "}; if (substr($1,1,1)!="#") print s;}' | sort
