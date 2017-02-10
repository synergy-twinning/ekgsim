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

