/*  
    Copyright (c) 2017 Institute Jožef Stefan, Jamova cesta 39, SI-1000, Ljubljana, Slovenija

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    Please cite the following works (bibtex source below):
        - DepolliAvbeljTrobec2008 for the simulator and simulation-based optimization
        - DepolliTrobecFilipic2013 for the AMS-DEMO optimizer
        - TrobecDepolliAvbelj2009 for the simulator

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
*/

#ifndef SIM_H_INCLUDED
#define SIM_H_INCLUDED


#include <AMS-DEMO/ParallelNumericOptimizer.h>
#include <AMS-DEMO/VectorArithmetics.h>
#include <AMS-DEMO/GeneralOptimizationAlgorithm.h>


/*******************************************************************************************//**
    class OutputSettings
    
    defines the data that should be exported as the result
**/
struct OutputSettings {
	std::vector<size_t> outputCellAps;
	bool layerAps;
	bool excitationSequence, result;
	
	OutputSettings() : layerAps(false), excitationSequence(false), result(false) {}
};


/*******************************************************************************************//**
    hidden class SimImplementation (defined in .cpp file)
**/
class SimImplementation;


/*******************************************************************************************//**
    implementation of VirtualOptimizationFunction
**/
struct OptimizationFunction : public VirtualOptimizationFunction {	
//	typedef std::vector<double> Value;
//	typedef std::vector<double> Input;
	SimImplementation* impl;
	size_t valueLen;
	size_t propertiesLen;
	std::vector<double> geneMin;
	std::vector<double> geneMax;

public:
    /// OptimizationFunction takes single parameter at initialization - the number of properties
	OptimizationFunction(size_t propertiesL);
	~OptimizationFunction();
	
	/// setup only sets the outputs
	void setup(const OutputSettings& outSet);
	/// operator() runs the optimization function on the given input (solution).
	/// Results are saved to parameters result, violation, and properties; function does not return anything
	/// if result and properties are vectors, they are automatically resized to correct length
	void operator() (const Input& solution, Value& result, double& violation, Properties& properties) const;
	/// query for parameters that define individual solutions and are located in implementation (number of genes, number of criteria, min and max for genes)	
	void getGeneParams(size_t& numGenes, size_t& numCriteria, std::vector<double>& gMin, std::vector<double>& gMax);
	/// query for parameters that define optimization and are located in implementation (number of generations, size of the population, length of queue)	
	void getEvolutionaryParams(size_t& numGenerations, size_t& popSize, int& queueSize);
	/// normalization of parameters
	void normalize(Input& solution) const;
};


void testWohlInterpolation();


#endif // SIM_H_INCLUDED
