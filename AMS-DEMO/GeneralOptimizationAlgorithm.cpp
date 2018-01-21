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

#include "GeneralOptimizationAlgorithm.h"
#include "VectorArithmetics.h"
#include "internalEval.h"
#include <stdexcept>


namespace {
	
	template<class T>
	void clip(T& variable, const T& minimum, const T& maximum) {
		if (variable > maximum)
			variable = maximum;
		else if (variable < minimum)
			variable = minimum;
	}
	
}


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// class Settings
///
Settings::Settings() {
    evaluation.internalFunctionDelay = 0;
	evaluation.inputFileName = "input.txt";
	evaluation.outputFileName = "output.txt";
	
	algorithmParams.deSchema = "rand/1/bin";
	algorithmParams.scalingFactors.resize(1, 0.5);
	algorithmParams.minQueueLength = algorithmParams.maxQueueLength = 1;
}


bool Settings::readIni(const char* fname) {
	Ini::File ini(fname);
	
	if (ini) {
		bool ok = 
			readEvaluationSection(ini) && 
			readOptimizationSection(ini) &&
			readInitialPopulationSection(ini);
		return ok;
	} else
		return false;
}


bool Settings::readEvaluationSection(Ini::File& ini) {
	int section = ini.getSectionNumber("evaluation");
	if (section != -1) {
	    int tempInt = 0;
		evaluation.debugOutput = (ini.loadVar(tempInt, "debug output", section) && (tempInt != 0));
		ini.loadVar(evaluation.internalFunctionName, "internal function", section);
		ini.loadVar(evaluation.internalFunctionDelay, "internal function delay", section);
		ini.loadVar(evaluation.inputFileName, "input file name", section);
		ini.loadVar(evaluation.outputFileName, "output file name", section);
		ini.loadVar(evaluation.commandLine, "command line", section);
		ini.loadVar(evaluation.chromosomeVectorLength, "chromosome vector length", section);
		ini.loadVar(evaluation.criteriaVectorLength, "criteria vector length", section);
		ini.loadVar(evaluation.propertiesVectorLength, "properties vector length", section);
		return true;
	} else {
		return false;
	}
}


bool Settings::readOptimizationSection(Ini::File& ini) {
	int section = ini.getSectionNumber("optimization");
	if (section != -1) {
		ini.loadVar(algorithmParams.randomSeed, "random seed", section);
		ini.loadVar(algorithmParams.populationSize, "population size", section);
		size_t numGenerations = 0;
		if (ini.loadVar(numGenerations, "max number of generations", section))
			algorithmParams.maxNumOfEvaluations = numGenerations * algorithmParams.populationSize;
		ini.loadVar(algorithmParams.maxNumOfEvaluations, "max number of evaluations", section);
		ini.loadVar(algorithmParams.deSchema, "DE schema", section);
		ini.loadVar(algorithmParams.crossoverProbability, "p crossover", section);
		{
			Ini::ArrayReader<double, std::vector<double> > scalingFactorsReader(algorithmParams.scalingFactors, 3);
			ini.loadVar(scalingFactorsReader, "scaling factors", section);
			if (algorithmParams.scalingFactors.size() == 0)
				algorithmParams.scalingFactors.push_back(0.5);
		}
		{
			std::string temp;
			if (ini.loadVar(temp, "queue length", section)) {
				std::istringstream tempStr(temp);
				tempStr >> algorithmParams.minQueueLength;
				algorithmParams.maxQueueLength = algorithmParams.minQueueLength;
				tempStr.ignore(1);
				if (tempStr) {
					tempStr >> algorithmParams.maxQueueLength;
				}
				if (algorithmParams.maxQueueLength < algorithmParams.minQueueLength)
					algorithmParams.maxQueueLength = algorithmParams.minQueueLength;
				clip(algorithmParams.minQueueLength, 1, 100);
				clip(algorithmParams.maxQueueLength, 1, 100);
			} else {
				algorithmParams.maxQueueLength = algorithmParams.minQueueLength = 1;
			}
		}
		return true;
	} else {
		return false;
	}
}


bool Settings::readInitialPopulationSection(Ini::File& ini) {
	int section = ini.getSectionNumber("initial population");
	if (section != -1) {
		ini.loadVar(initialPopulation.readFilename, "file", section);
		initialPopulation.readGeneration = -1;
		ini.loadVar(initialPopulation.readGeneration, "read generation", section);
		initialPopulation.randomRatio = 0;
		ini.loadVar(initialPopulation.randomRatio, "min random ratio", section);
		
		{
			Ini::ArrayReader<double, std::vector<double> > reader(initialPopulation.geneMin);
			ini.loadVar(reader, "gene min", section);
			if (initialPopulation.geneMin.size() == 0)
				initialPopulation.geneMin.push_back(0);
		}
		{
			Ini::ArrayReader<double, std::vector<double> > reader(initialPopulation.geneMax);
			ini.loadVar(reader, "gene max", section);
			if (initialPopulation.geneMax.size() == 0)
				initialPopulation.geneMax.push_back(1);
		}
		return true;
	} else {
		return false;
	}
}


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// class GeneralOptimizationFunction
///
GeneralOptimizationFunction::Properties GeneralOptimizationFunction::dummyProperties;


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// class GeneralOptimizationAlgorithm
///
bool GeneralOptimizationAlgorithm::loadSettings(const char* fname) {
	if (settings.readIni(fname)) {
		return updateSettings();
	} else {
		// ini could not be read
		return false;
	}
}


bool GeneralOptimizationAlgorithm::updateSettings() {			
	if (settings.evaluation.criteriaVectorLength == 1)
		eval = &soAlg.evaluation();
	else
		eval = &moAlg.evaluation();
	
	if (settings.evaluation.internalFunctionName != "") {
		// internal evaluation - choose from test functions implemented in pDEMO or supplied as shared libraries;
		// name must be two part (separated with space), first part declares function set (like CEC2007), and secon part
		// the function name within that set (like OKA2)
		std::istringstream ostr(settings.evaluation.internalFunctionName);
		std::string part;
		std::getline(ostr, part, ' ');
		if (settings.evaluation.debugOutput)
            std::cerr << "internal function selected from family " << part;
		if (part == "CEC2007") {
            std::shared_ptr<InternalOptimizationFunction<TestFunctionCEC2007> > virtPtr(
                new InternalOptimizationFunction<TestFunctionCEC2007> (
                    settings.evaluation.criteriaVectorLength,
                    settings.evaluation.propertiesVectorLength,
                    settings.evaluation.debugOutput
                    ));
            std::getline(ostr, part, ' ');
            if (settings.evaluation.debugOutput)
                std::cerr << ", function \n" << part;
            virtPtr->chooseFuncByName(part.c_str(), settings.evaluation.internalFunctionDelay);
            settings.evaluation.internalFunction = virtPtr;
		} else {
		    if (settings.evaluation.debugOutput)
                std::cerr << "\n";
        }
            
		if (settings.evaluation.internalFunction.get() != 0) {
			eval->setImplementation(settings.evaluation.internalFunction);
		} else
			throw std::runtime_error("missing implementation for internal evaluation function");
	} else {
		ExternalOptimizationFunction* tempPtr = new 
			ExternalOptimizationFunction(	settings.evaluation.criteriaVectorLength,
											settings.evaluation.propertiesVectorLength,
                                            settings.evaluation.debugOutput);
		
		tempPtr->setNormalization(	settings.initialPopulation.geneMin,
									settings.initialPopulation.geneMax);
		
#ifdef NO_MPI
		tempPtr->setCommand(
			settings.evaluation.commandLine, 
			"",
			settings.evaluation.inputFileName,
			settings.evaluation.outputFileName);
#else
		std::ostringstream tempSS;
		if (communicator.getSize() != 1)
			tempSS << "process" << communicator.getRank();
		tempPtr->setCommand(
			settings.evaluation.commandLine, 
			tempSS.str().c_str(),
			settings.evaluation.inputFileName,
			settings.evaluation.outputFileName);
#endif //NO_MPI
		
		eval->setImplementation(tempPtr);
	} 
	
	soAlg.debugOutput = settings.evaluation.debugOutput;
	soAlg.randomSeed = settings.algorithmParams.randomSeed;
	soAlg.scalingFactor.resize(settings.algorithmParams.scalingFactors.size());
	std::copy(settings.algorithmParams.scalingFactors.begin(), 
		settings.algorithmParams.scalingFactors.end(),
		soAlg.scalingFactor.begin());
	soAlg.crossoverProbability = settings.algorithmParams.crossoverProbability;
	soAlg.maxNumOfEvaluations = settings.algorithmParams.maxNumOfEvaluations;
	soAlg.selectSchema(settings.algorithmParams.deSchema);
	
	moAlg.debugOutput = settings.evaluation.debugOutput;
	moAlg.randomSeed = settings.algorithmParams.randomSeed;
	moAlg.scalingFactor.resize(settings.algorithmParams.scalingFactors.size());
	std::copy(settings.algorithmParams.scalingFactors.begin(), 
		settings.algorithmParams.scalingFactors.end(),
		moAlg.scalingFactor.begin());
	moAlg.crossoverProbability = settings.algorithmParams.crossoverProbability;
	moAlg.maxNumOfEvaluations = settings.algorithmParams.maxNumOfEvaluations;
	moAlg.selectSchema(settings.algorithmParams.deSchema);
	
	return true;
}


void GeneralOptimizationAlgorithm::run() {
	// generation of initial population
	{	
		DynamicRandomInitializer<double, Individual> randomInit(
			settings.evaluation.chromosomeVectorLength, settings.evaluation.criteriaVectorLength);
		
		if (settings.initialPopulation.readFilename != "") {
			try {
				IndividualsFile file(settings.initialPopulation.readFilename.c_str());
				file.readGeneration(randomInit.preLoaded, settings.initialPopulation.readGeneration);
			} catch (const char*) {
				// ignore custom exceptions which only mean the file could not be read
			} catch (...) {
				throw;
			}
		}
		
		for (int i = randomInit.preLoaded.size() - (int)ceil(
			(1.0f - settings.initialPopulation.randomRatio) * settings.algorithmParams.populationSize);
			i > 0; --i) 
		{
			size_t d = Random::CRand().exclusiveInterval(randomInit.preLoaded.size());
			if ((d+1) < randomInit.preLoaded.size()) 
				randomInit.preLoaded[d] = randomInit.preLoaded.back();
		
			randomInit.preLoaded.pop_back();
		}
		
		for (size_t i = 0; i < settings.evaluation.chromosomeVectorLength; ++i) {
			randomInit.defineGeneValueRange(i, 
				settings.initialPopulation.geneMin[std::min(i, settings.initialPopulation.geneMin.size()-1)], 
				settings.initialPopulation.geneMax[std::min(i, settings.initialPopulation.geneMax.size()-1)]);
		}
		if (settings.evaluation.criteriaVectorLength == 1) {
			soAlg.initPopulation(settings.algorithmParams.populationSize, randomInit);
			soAlg.minQueueLength = settings.algorithmParams.minQueueLength;
			soAlg.maxQueueLength = settings.algorithmParams.maxQueueLength;
		} else {
			moAlg.initPopulation(settings.algorithmParams.populationSize, randomInit);
			moAlg.minQueueLength = settings.algorithmParams.minQueueLength;
			moAlg.maxQueueLength = settings.algorithmParams.maxQueueLength;
		}
	}
	
	if (settings.evaluation.criteriaVectorLength == 1) {
		// single objective
#ifdef NO_MPI
		soAlg.evolution();
#else
		soAlg.evolution(communicator);
#endif //NO_MPI
	} else {
		// multi objective
#ifdef NO_MPI
		moAlg.evolution();
#else
		moAlg.evolution(communicator);
#endif //NO_MPI
	}
}

#ifndef NO_MPI	
void GeneralOptimizationAlgorithm::runWorker() {
	settings.evaluation.criteriaVectorLength == 1 ?
		soAlg.evolution(communicator) :
		moAlg.evolution(communicator);
}
#endif //NO_MPI


const std::vector<GeneralOptimizationAlgorithm::Individual>& GeneralOptimizationAlgorithm::population() const {
	return (settings.evaluation.criteriaVectorLength == 1 ? soAlg.population() : moAlg.population());
}


void GeneralOptimizationAlgorithm::getFront(std::vector<size_t>& front) const {
	if (settings.evaluation.criteriaVectorLength == 1)
		throw std::runtime_error("Programming stupidity detected in [GeneralOptimizationAlgorithm::getFront]: single objective algorithm does not have a front");
	moAlg.getFront(front);
}
