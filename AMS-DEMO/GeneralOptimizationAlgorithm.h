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

#ifndef GENERALOPTIMIZATIONALGORITHM_H_INCLUDED
#define GENERALOPTIMIZATIONALGORITHM_H_INCLUDED

#ifndef NO_MPI
//#include "jobDistributer.h"
#endif

#include "ParallelNumericOptimizer.h"
#include "pDe.h"
#include "pDemo.h"
#include "TypeWrapper.h"

#include <Ini.h>
#include <vector>
#include <string>
#include <memory>


struct VirtualOptimizationFunction {
	typedef std::vector<double> Value;
	typedef std::vector<double> Input;
	typedef std::vector<double> Properties;
	bool geneBounds;
	
	VirtualOptimizationFunction() : geneBounds(true) {}
	virtual ~VirtualOptimizationFunction() {}
	virtual void operator() (const Input& solution, Value& result, double& violation, Properties& properties) const = 0;
	virtual void normalize(Input& solution) const {}
};


class Settings {
public:
	struct {
	    bool debugOutput;
		std::string internalFunctionName;
		int internalFunctionDelay;
		std::shared_ptr<VirtualOptimizationFunction> internalFunction;
		std::string inputFileName, outputFileName, commandLine;
		size_t chromosomeVectorLength, criteriaVectorLength, propertiesVectorLength;
	} evaluation;
	
	struct {
		int randomSeed;
		size_t populationSize;
		size_t maxNumOfEvaluations;
		double crossoverProbability;
		std::vector<double> scalingFactors;
		std::string deSchema;
		// parallel algorithm only:
		int minQueueLength, maxQueueLength;
	} algorithmParams;
	
	struct {
		std::vector<double> geneMin;
		std::vector<double> geneMax;
		float randomRatio;
		int readGeneration;
		std::string readFilename;
	} initialPopulation;
	
public:
	Settings();
	bool readIni(const char* fname);
	
protected:
	bool readEvaluationSection(Ini::File& ini);
	bool readOptimizationSection(Ini::File& ini);
	bool readInitialPopulationSection(Ini::File& ini);
};


struct Normalization {
    std::vector<double> geneMin;
    std::vector<double> geneMax;
    
    void normalize(std::vector<double>& in) const {
        for (size_t i = 0; (i < geneMin.size()) && (i < in.size()); ++i)
            if (in[i] < geneMin[i]) in[i] = geneMin[i];
        for (size_t i = 0; (i < geneMax.size()) && (i < in.size()); ++i)
            if (in[i] > geneMax[i]) in[i] = geneMax[i];
    }
};


/// wrapper of ExternalEvaluation template class
struct ExternalOptimizationFunction : public VirtualOptimizationFunction {
	ExternalEvaluation<Input, Value, Properties, Normalization> eval;
	size_t valueLen;
	size_t propertiesLen;
	bool debugEnabled;
	
public:
	ExternalOptimizationFunction(size_t valueL, size_t propertiesL, bool debug = true) {
		valueLen = valueL;
		propertiesLen = propertiesL;
		debugEnabled = debug;
	}
	
	void setCommand(const std::string& com, const std::string& dir = "", const std::string& inFile = "input.txt", const std::string& outFile = "output.txt") {
		eval.executableCommand = com;
		eval.inFname = inFile;
		eval.outFname = outFile;
		eval.homeDir = dir;
	}
	
	void setNormalization(const std::vector<double>& geneMin, const std::vector<double>& geneMax) {
		eval.geneMin = geneMin;
		eval.geneMax = geneMax;
	}
	
	void operator() (const Input& solution, Value& result, double& violation, Properties& properties) const {
		static size_t num = 0;
		++num;
		if (debugEnabled)
            std::cerr << "eval #" << num << "\r" << std::flush;
		if (valueLen > 0) TypeWrapper::resize(result, valueLen);
		if (propertiesLen > 0) TypeWrapper::resize(properties, propertiesLen);
		eval(solution, result, violation, properties);
	}
	
	void normalize(Input& solution) const {
		eval.normalize(solution);
	}
};


/// wrapper of InternalEvaluation template class
template<class FuncClass>
struct InternalOptimizationFunction : public VirtualOptimizationFunction {
	size_t valueLen;
	size_t propertiesLen;
	FuncClass func;
	Normalization normalization;
	bool debugEnabled;
	
public:
	InternalOptimizationFunction(size_t valueL, size_t propertiesL = 0, bool debug = true) {
		valueLen = valueL;
		propertiesLen = propertiesL;
		debugEnabled = debug;
	}
	
	void setNormalization(const std::vector<double>& geneMin, const std::vector<double>& geneMax) {
		normalization.geneMin = geneMin;
		normalization.geneMax = geneMax;
	}
	
	void chooseFuncByName(const char* name, int delay = 0) {
	    func.selectFunction(name);
	    func.setSleepTime(delay);
	}
	
	void operator() (const Input& solution, Value& result, double& violation, Properties& properties) const {
		static size_t num = 0;
		++num;
		if (debugEnabled) {
            std::cerr << "eval #" << num << "\r";
            if (num && 0x0F == 0)
                std::cerr << std::flush;
            }
		if (valueLen > 0) TypeWrapper::resize(result, valueLen);
		if (propertiesLen > 0) TypeWrapper::resize(properties, propertiesLen);
		func(solution, result, violation, properties);
	}
	
	void normalize(Input& solution) const {
		normalization.normalize(solution);
	}
};


struct GeneralOptimizationFunction {
	typedef std::vector<double> Value;
	typedef std::vector<double> Input;
	typedef std::vector<double> Properties;
	static Properties dummyProperties;
	
	bool geneBounds;
	std::shared_ptr<VirtualOptimizationFunction> functionImplementation;
	
	// implementarion must be allocated with operator new!
	void setImplementation(VirtualOptimizationFunction* imp) {
		functionImplementation = std::shared_ptr<VirtualOptimizationFunction>(imp);
		geneBounds = functionImplementation->geneBounds;
	}
	
	void setImplementation(std::shared_ptr<VirtualOptimizationFunction>& imp) {
		functionImplementation = imp;
		geneBounds = functionImplementation->geneBounds;
	}
	
	inline void operator() (const Input& solution, Value& result, double& violation, Properties& properties = dummyProperties) const {
		assert(&*functionImplementation != 0);
		(*functionImplementation)(solution, result, violation, properties);
	}
	
	void normalize(Input& solution) const {
		assert(&*functionImplementation != 0);
		functionImplementation->normalize(solution);
	}
};


class GeneralOptimizationAlgorithm {
public:
	typedef std::vector<double> 			Chromosome;
	typedef PDemo<Chromosome, GeneralOptimizationFunction> 	MultiObjectiveAlg;
	typedef PDe<Chromosome, GeneralOptimizationFunction> 	SingleObjectiveAlg;
	typedef SingleObjectiveAlg::Individual	Individual;
//	typedef IndividualGenerator<GeneralOptimizationFunction::Input,
//		GeneralOptimizationFunction::Value>					IndividualGenerator;
	
public:
	Settings settings;
	
private:
	SingleObjectiveAlg 	soAlg;
	MultiObjectiveAlg	moAlg;
	GeneralOptimizationFunction* eval;
	
#ifndef NO_MPI
	Mpi::Communicator communicator;
#endif //NO_MPI

public:
	// load settings from file and call updateSettings
	bool loadSettings(const char* fname);
	// when settings are set, call this function to update algorithm accordingly
	// (settings are just a common and user friendly placeholder, algorithm must
	// be set explicitly)
	bool updateSettings();
	
#ifndef NO_MPI
	void setCommunicator(const Mpi::Communicator& com) {
		communicator = com;
	}
	
	void runWorker();
#endif //NO_MPI
	
	void run();
	const std::vector<Individual>& population() const;
	void getFront(std::vector<size_t>& front) const;
};


#endif // GENERALOPTIMIZATIONALGORITHM_H_INCLUDED

