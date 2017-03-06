#!/bin/bash

# takes individuals file as an input and outputs hypervolume on std.out
echo ${1-"individuals.txt"} 
../AMS-DEMO/DEMO -individuals:${1-"individuals.txt"} -analysis:hypervolume -gen:0:-1 > hyper.txt
cat hyper.txt
