# README #

This repository contains the source for ECG simulator based on 3D voxel model of a human heart, analytic model of action potentials and cellular automaton for excitation sequence. It contains embedded optimization for determining AP model parameters from the prescribed ECG shape, the given heart model, and the given ECG measuring points (lead positions)

# Set-up

## Preparing the environment

Several operating system packages are required to compile and run EkgSim. 
These are usually already provided on HPC but might are not installed by default on all Linux systems. 
An example of how the environment should be prepared on a new Ubuntu installation is given below:

~~~
apt-get update
apt-get install make build-essential libopenmpi-dev unzip openmpi-bin ssh
~~~

If you want to download the latest stable version (original simulator) use:
~~~
wget https://github.com/synergy-twinning/ekgsim/archive/v0.1.5.zip -O ekgsim.zip
unzip ekgsim.zip
cd ekgsim-0.1.5/
~~~

If you want to download the latest master version from the repository use:
~~~
wget https://github.com/synergy-twinning/ekgsim/archive/master.zip -O ekgsim.zip
unzip ekgsim.zip
cd ekgsim-master/
~~~

## Compiling the code

Make will create an executable in the directory Release. This is an MPI-enabled parallel program but can be run sequentially if needed.

Run `make -f makefile.manual` inside the ekgSim directory.
Note that makefiles were generated by CodeLite in Ubuntu environment and might contain several absolute paths. These should in principle not matter when running make. If make fails, however, see if it is due to absolute/relative paths or something being miss-configured in `ekgSim.mk`, `simlib/simlib.mk` or `copyOfLibs/copyOfLibs.mk`. Currently we do not have the knowledge or resources (e.g. time) to enable magic build on all the platforms, so a bit of configuration is to be expected.

You will probably see some warnings during the compilation. Don't panic, everything should still work at the end.
```
In file included from Ini.cpp:65:0:
Ini.h: In constructor 'Ini::File::File(const char*)':
Ini.h:146:36: warning: 'Ini::File::separators_' will be initialized after [-Wreorder]
         std::string                separators_;
                                    ^
Ini.h:144:36: warning:   'std::string Ini::File::filename_' [-Wreorder]
         std::string                filename_;
                                    ^
Ini.cpp:209:2: warning:   when initialized here [-Wreorder]
  File::File(const char* filename)
  ^
Ini.cpp: In member function 'void Ini::File::parse(const string&, std::string&, std::string&) const':
Ini.cpp:368:26: warning: comparison between signed and unsigned integer expressions [-Wsign-compare]
                 if (i == firstChar)
                          ^
```

# Running the simulator

After a successful make, the executable is generated as *Release/ekgSim*.
The simulator requires an environment to be set up. An example of such environment is given in directory *testRun*. The simplest way to running the executable is to copy the executable to this directory and then run everything from there:
~~~
cp Release/ekgSim testRun/
cd testRun
~~~

## Sequential execution

### Versions < v0.1.4?
The -sim parameter is a vector of numeric parameters of length 21, separated by commas (no spaces!). Parameters need to be in the following ranges:
~~~
min = 1.5, 0.85, 0.05, 0.0003, 0.01,  0.01, 200, 1.5, 0.85, 0.05, 0.0003, 0.01,  0.01, 200, 1.5, 0.85, 0.05, 0.0003, 0.01,  0.01, 200
max = 3.5, 0.95, 0.2,  0.0010, 0.10,  0.10, 400, 3.5, 0.95, 0.2,  0.0010, 0.10,  0.10, 400, 3.5, 0.95, 0.2,  0.0010, 0.10,  0.10, 400
~~~

Here is one example on how to call the simulator sequentially.

~~~~
./ekgSim test -sim 2.25214,0.925452,0.0943991,0.00035813,0.0890636,0.0632915,226.183,2.08002,0.857738,0.162565,0.000369406,0.0965625,0.0523254,232.278,3.07406,0.855509,0.125351,0.000710767,0.0720323,0.0187579,200.93 -out result
~~~~
The result of the simulation in form of a simulated ECG signal will be stored in _result.column_ file, which is in Matlab/octave readable text format.
File _result.column_ will contain 3 columns, first is the time in milliseconds, followed by the 2 simulated ECGs for all the measuring points.

Beside the output file *result*, some of the results are also output to the standard output, for example:
~~~~
...
eval 1   simulation done in 439.131 seconds
 criteria = <0.099919,0.012262>, violation = 0
All done
~~~~

In order to use the simulator for the optimization, you will need to parse the standard output for the criteria.

### Versions >= v0.1.4
In this version the simulator was modified to accept additional parameters for electrodes positioning. Parameters now include:

- 4 x 3 (5th, 6th, 7th and 8th) wohlfart's parameters
- 2 x 2 electrodes position parameters

Together that is 16 parameters with the following bounds:
 ~~~
min = 0.0003, 0.01,  0.01, 200, 0.0003, 0.01,  0.01, 200, 0.0003, 0.01,  0.01, 200, -50, -50, -50, -50
max = 0.0010, 0.10,  0.10, 400, 0.0010, 0.10,  0.10, 400, 0.0010, 0.10,  0.10, 400, 50, 50, 50, 50
~~~

Use the same call as before to run the simulation, for example:
~~~~
./ekgSim test -sim 0.00035813,0.0890636,0.0632915,226.183,0.000369406,0.0965625,0.0523254,232.278,0.000710767,0.0720323,0.0187579,200.93,23,22,15,13 -out result
~~~~

# Running the optimization
## Sequential mode

~~~~
# execute the simulator without extra parameters
cp Release/ekgSim testRun/
cd testRun
./ekgSim
~~~~

Several files will be created as a result: 
- _individuals.txt_, which contains individuals that reside in the working population of the evolutionary optimization algorithm 
- _evaluations.txt_, which contains all the performed evaluations (even the ones that result in individuals with high (high = bad) value of objective functions)
- _featureFile.txt_, which contains the details of the execution (host name(s), execution time, etc)

## Parallel mode

Parallel execution is supported on MPI enabled systems. Parallelization is done on the level of optimization only. This means that optimization will evaluate individuals (run individual simulations and evaluate their results under one or more criteria) in parallel on multiple machines or processes within the same machine. Apart from the execution command line, everything is the same as when running the program sequentially.

~~~~
# execute the simulator through MPI call on 16 processors:
cp Release/ekgSim testRun/
cd testRun
mpirun -n 16 ekgSim
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
