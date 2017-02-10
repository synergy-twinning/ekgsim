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
	// not yet used:
	struct Normalization {
		std::vector<double> geneMin;
		std::vector<double> geneMax;
		
		void normalize(Input& in) const {
			for (size_t i = 0; (i < geneMin.size()) && (i < in.size()); ++i)
				if (in[i] < geneMin[i]) in[i] = geneMin[i];
			for (size_t i = 0; (i < geneMax.size()) && (i < in.size()); ++i)
				if (in[i] > geneMax[i]) in[i] = geneMax[i];
		}
	};
	
//	typedef std::vector<double> Value;
//	typedef std::vector<double> Input;
	SimImplementation* impl;
	size_t valueLen;
	size_t propertiesLen;

public:
    /// PptimizationFunction takes single parameter at initialization - the number of properties
	OptimizationFunction(size_t propertiesL);
	~OptimizationFunction();
	
	/// setup only sets the outputs
	void setup(const OutputSettings& outSet);
	/// operator() runs the optimization function on the given input (solution).
	/// Results are saved in to result, violation, and properties; function does not return anything
	/// if result and properties are vectors, they are automatically resized to correct length
	void operator() (const Input& solution, Value& result, double& violation, Properties& properties) const;
	/// query for parameters that define individual solutions and are located in implementation (number of genes, number of criteria, min and max for genes)	
	void getGeneParams(size_t& numGenes, size_t& numCriteria, std::vector<double>& gMin, std::vector<double>& gMax);
	/// query for parameters that define optimization and are located in implementation (number of generations, size of the population, length of queue)	
	void getEvolutionaryParams(size_t& numGenerations, size_t& popSize, int& queueSize);
};


void testWohlInterpolation();


#endif // SIM_H_INCLUDED
