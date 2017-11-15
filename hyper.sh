#!/bin/bash

# takes an optional argument - the name of input file
inFile=${1-"testRun/individuals.txt"}

# takes individuals file as an input and outputs hypervolume on std.out
echo -e "using input file:\n   $inFile"
echo -e "output will be written to file:\n   hyper.txt"

./AMS-DEMO/DEMO -individuals:${inFile} -analysis:hypervolume -params:0,2:0,2 -gen:0:-1 > hyper.txt
# optionally, you might wish to let AMS-DEMO calculate the limits for hypervolume by itself; this can be done by using (uncomment it) the line below:
#./AMS-DEMO/DEMO -individuals:${inFile} -analysis:hypervolume -gen:0:-1 > hyper.txt

echo ""
