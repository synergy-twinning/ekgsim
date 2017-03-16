# README #

This repository contains the source for ECG simulator based on 3D voxel model of a human heart, analytic model of action potentials and cellular automaton for excitation sequence. It contains embedded optimization for determining AP model parameters from the prescribed ECG shape, the given heart model, and the given ECG measuring points (lead positions)

# Compiling the code

Run `make -f makefile.manual` inside the ekgSim directory.
Note that makefiles was generated by CodeLite in Ubuntu environment and might contain several absolute paths. These should in principle not matter when running make. If make fails, however, see if it is due to absolute paths or something being miss-configured in ekgSim.mk, simlib/simlib.mk, or copyOfLibs/copyOfLibs.mk

## Preparing the environment

Several operating system packages are required to compile and run EkgSim. 
These are usually already provided on HPC but might are not installed by default on all Linux systems. 
An example of how the environment should be prepared on a new Ubuntu installation is given below:

~~~
apt-get update
apt-get install make build-essential libopenmpi-dev unzip openmpi-bin ssh
wget https://github.com/synergy-twinning/ekgsim/archive/master.zip -O ekgsim.zip
unzip ekgsim.zip
cd ekgsim-master
~~~

Make will create an executable in directory Release. This is an MPI-enabled parallel program but can be run sequentially if needed.

# Running the simulator

After a successful make, the executable is generated as *Release/ekgSim*.
The simulator requires an environment to be set up. An example of such environment is given in directory *testRun*. The simplest way to running the executable is to copy the executable to this directory and then run everything from there:
~~~
cp Release/ekgSim testRun/
cd testRun
~~~
Before running the executable, adjust the parameters of simulator and optimization in *simulator.ini*.
Also set up the linked files (filenames are specified in simulator.ini).

## Sequential execution

~~~~
# parameters is a vector of numeric parameters (of the correct length, which is set in simulator.ini), separated by commas (no spaces!)
Release/ekgSim test -sim parameters -out result
# the result of the simulation will be stored in result.column, which is in Matlab/octave readable text format
# result.column will contain _n+1_ columns, first is the time in milliseconds, followed by the _n_ simulated ECGs for all the measuring points
~~~~

Some of the results are also output to the standard output, for example:
~~~~
...
eval 1   simulation done in 439.131 seconds
 criteria = <0.099919,0.012262>, violation = 0
All done
~~~~

# Running the optimization
## Sequential mode

~~~~
# execute the simulator without extra parameters
Release/ekgSim
~~~~

Several files will be created as a result: 
- individuals.txt, which contains individuals that reside in the working population of the evolutionary optimization algorithm 
- evaluations.txt, which contains all the performed evaluations (even the ones that result in individuals with high (high = bad) value of objective functions)
- featureFile.txt, which contains the details of the execution (host name(s), execution time, etc)

## Parallel optimization

Parallel execution is supported on MPI enabled systems. Parallelization is done on the level of optimization only. This means that optimization will evaluate individuals (run individual simulations and evaluate their results under one or more criteria) in parallel on multiple machines or processes within the same machine. Apart from the execution command line, everything is the same as when running the program sequentially.

~~~~
# execute the simulator through MPI call on 16 processors:
mpirun -n 16 Release/ekgSim
~~~~

# Citation
Please cite the following works (bibtex source below):

- _DepolliAvbeljTrobec2008_ for the simulator and simulation-based optimization
- _DepolliTrobecFilipic2013_ for the AMS-DEMO optimizer
- _TrobecDepolliAvbelj2009_ for the simulator

```
@article{DepolliAvbeljTrobec2008,
    author = {Depolli, Matjaž and Avbelj, Viktor and Trobec, Roman},
    title = {Computer-Simulated Alternative Modes of {U}-Wave Genesis},
    journal = {Journal of Cardiovascular Electrophysiology},
    volume = {19},
    number = {1},
    publisher = {Blackwell Publishing Inc},
    issn = {1540-8167},
    url = {http://dx.doi.org/10.1111/j.1540-8167.2007.00978.x},
    doi = {10.1111/j.1540-8167.2007.00978.x},
    pages = {84--89},
    keywords = {U wave, ECG, action potential, repolarization, myocardium, computer simulation},
    year = {2008}
}

@article{DepolliTrobecFilipic2013,
    author = {Depolli, Matjaž and Trobec, Roman and Filipič, Bogdan},
    title = {Asynchronous master-slave parallelization of differential evolution for multiobjective optimization},
    journal = {Evolutionary Computation},
    volume = {21},
    number = {2},
    pages = {261-291},
    doi = {10.1162/EVCO_a_00076},
    issn = {1063-6560},
    url = {http://www.mitpressjournals.org/doi/abs/10.1162/EVCO_a_00076},
    year = {2013}
}

@inproceedings{TrobecDepolliAvbelj2009,
    title           = {Simulation of {ECG} repolarization phase with improved model of cell action potentials},
    author          = {Trobec, Roman and Depolli, Matja{\v{z}} and Avbelj, Viktor},
    booktitle       = {International Joint Conference on Biomedical Engineering Systems and Technologies},
    pages           = {325--332},
    year            = {2009},
    organization    = {Springer}
}
```
