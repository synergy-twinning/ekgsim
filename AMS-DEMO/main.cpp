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

#include "pch.h"

/*
 * everything below was moved to pch.h (a precompiled header)
 * 
#ifndef NO_MPI
	#include "mpiWrapper.h"
#endif
#include "VectorArithmetics.h"
#include "ParallelNumericOptimizer.h"
#include "pDemo.h"
#include "TestFunctions.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>
*/
#include "VectorArithmetics.h"
#include "ParallelNumericOptimizer.h"
#include "pDemo.h"
#include "TestFunctions.h"
#include <Arguments.h>
#include "GeneralOptimizationAlgorithm.h"
#include "HyperVolume.h"



namespace UnitTests {
	
	void testIteratorRandomizer() {
		std::vector<size_t> indices;
		randomizeIndices(indices, 100);
		std::cout << "\n";
		std::copy(indices.begin(), indices.end(), std::ostream_iterator<int> (std::cout, " "));
		std::cout << "\n";
		std::sort(indices.begin(), indices.end());
		std::copy(indices.begin(), indices.end(), std::ostream_iterator<int> (std::cout, " "));
		std::cout << "\n";
	}

	void testPDemo() {
		static const int dimensionality = 30;
		typedef TestFunctions::ZDT1<dimensionality> Problem;
		typedef Array<double, dimensionality> Chromosome;
		typedef PDemo<Chromosome, Problem> DemoAlg;
		DemoAlg alg;
		
		{
			DynamicRandomInitializer<double, DemoAlg::Individual> randomInit(dimensionality, 2);
			for (size_t i = 0; i < Problem::dimensionality; ++i)
				randomInit.defineGeneValueRange(i, 0.0, 1.0);
			alg.initPopulation(100, randomInit);
		}
		alg.scalingFactor.resize(2, 0.5f);
		alg.crossoverProbability = 0.3;
		alg.maxNumOfEvaluations = 10000;
		alg.evaluationsPerTruncation = 150;
		alg.evolution();
		
		
#ifndef NO_MPI
		if (Mpi::Communicator().getRank() == 0)
#endif //NO_MPI
		{
			{
				std::ofstream out("pop.txt");
				for (size_t i = 0; i < alg.population().size(); ++i) {
					out << alg.population()[i].criteria[0] << "\t" << alg.population()[i].criteria[1] 
						<< "\t" << i << "\t" << alg.population()[i].violation << "\t";
					
					for (size_t j = 0; j < alg.population()[i].chromosome.size(); ++j)
						out << alg.population()[i].chromosome[j] << "\t";
					out << "\n";
				}
			}
			{
				std::ofstream out("front.txt");
				std::vector<size_t> front;
				alg.getFront(front);
				for (size_t i = 0; i < front.size(); ++i) {
					out << alg.population()[front[i]].criteria[0] << "\t" << alg.population()[front[i]].criteria[1] << "\t";
					out << front[i] << "\t";
					for (size_t j = 0; j < alg.population()[front[i]].chromosome.size(); ++j)
						out << alg.population()[front[i]].chromosome[j] << "\t";
					out << "\n";
				}
			}
		}
	}
}


/// ************************************************************************************
/// optimization 
/// ************************************************************************************
#ifdef NO_MPI
void runOptimizer(const std::string& settingsFname)
#else
void runOptimizer(const std::string& settingsFname, Mpi::Environment& mpi)
#endif //NO_MPI
{
	GeneralOptimizationAlgorithm alg;
	if (alg.loadSettings(settingsFname.c_str())) {
#ifdef NO_MPI
		std::cout << "AMS DEMO v.3\n" << std::endl;
		alg.run();
#else	
		alg.setCommunicator(mpi.getCommunicator());
		
		if (mpi.getCommunicator().getRank() == 0) {
			std::cout << "AMS DEMO v.3 (" << mpi.getCommunicator().getSize() << " CPUs)\n" << std::endl;
			alg.run();
		} else {
			alg.runWorker();
		}
		if (mpi.getCommunicator().getRank() == 0) 
#endif
		{
			if (alg.settings.evaluation.criteriaVectorLength == 1) {		
				// in single objective version output best individual
				size_t bestI = 0;
				double bestVal = alg.population()[0].criteria[0];
				
				for (size_t i = 1; i < alg.population().size(); ++i) {
					if (alg.population()[i].criteria[0] < bestVal) {
						bestVal = alg.population()[i].criteria[0];
						bestI = i;
					}
				}
				std::cout << "minimum value = " << bestVal << ", " 
					<< alg.population()[bestI].getProperties() << " for " 
					<< alg.population()[bestI].chromosome << "\n";
			} else {
				// in multiobjective version output front as front.txt
				std::ofstream out("front.txt");
				std::vector<size_t> front;
				alg.getFront(front);
				for (size_t i = 0; i < front.size(); ++i) {
					out << alg.population()[front[i]].criteria[0];
					for (size_t j = 1; j < alg.population()[0].criteria.size(); ++j)
						out << '\t' << alg.population()[front[i]].criteria[j];
					out << '\n';
				}
				std::cout << "last front saved as front.txt\n";
			}
		}
	} else {
		std::cout << "could not open " << settingsFname << "\n";
	}
}


/// ************************************************************************************
/// analysis
/// ************************************************************************************
enum AnalysisMode {
	an_default = 0,
	an_front,
	an_hypervolume,
	an_best
};


enum OutputMode {
	out_all = 0,
	out_chromosome,
	out_violation,
	out_properties,
	out_criteria,
	out_rank,
	out_evaluation_time,
	out_life_time,
	out_no_gen
};


typedef IndividualsFile::FileLine<GeneralOptimizationAlgorithm::Individual> IndFileLine;


void outputIndividual(const IndFileLine& fline, size_t genNum, int o) {
    bool outGen = true;
    if ((out_no_gen & o) > 0) {
        outGen = false;
        o -= out_no_gen;
    }
    
    if (outGen) {
        std::cout << genNum << " \t" ;
    } else {
    }
    
	switch(o) {
	default:
	case out_all:
		std::cout << fline.individual.chromosome << " \t" 
			<< fline.individual.violation << " \t" << fline.individual.properties << " \t" 
			<< fline.individual.criteria;
		if (fline.rank != -1)
			std::cout << " \t" << fline.rank << " \t" << fline.evalTime << " \t" 	
				<< fline.lifetime;
		std::cout << "\n";
		break;
	case out_chromosome:
		std::cout << fline.individual.chromosome[0];
		for (size_t i = 1; i < fline.individual.chromosome.size(); ++i)
			std::cout << " " << fline.individual.chromosome[i];
		std::cout << "\n";
		break;
	case out_criteria:
		std::cout << fline.individual.criteria[0];
		for (size_t i = 1; i < fline.individual.criteria.size(); ++i)
			std::cout << " " << fline.individual.criteria[i];
		std::cout << "\n";
		break;
	case out_properties:
		if (GeneralOptimizationAlgorithm::Individual::hasProperties) {
			if (fline.individual.getProperties().size() > 0)
				std::cout << fline.individual.getProperties()[0];
			else
				std::cout << "/";
			for (size_t i = 1; i < TypeWrapper::size(fline.individual.getProperties()); ++i)
				std::cout << " " << fline.individual.getProperties()[i];
		}
		std::cout << "\n";
		break;
	case out_violation:
		std::cout << fline.individual.violation << "\n";
		break;
	case out_rank:
		std::cout << fline.rank << "\n";
		break;
	case out_evaluation_time:
		std::cout << fline.evalTime << "\n";
		break;
	case out_life_time:
		std::cout << fline.lifetime << "\n";
		break;
	}
}


void outputIndividual(const GeneralOptimizationAlgorithm::Individual& individual, size_t genNum, int o) {
	IndFileLine fline;
	fline.individual = individual;
	outputIndividual(fline, genNum, o);
}


struct OutputInd {
	int outMode;
	
	OutputInd(int m) : outMode(m) {}
	
	void operator() (const IndFileLine& fline, size_t genNum) {
		outputIndividual(fline, genNum, outMode);
	}
};


struct FormFront {
	int outMode;
	std::map<size_t, std::vector<GeneralOptimizationAlgorithm::Individual> > generations; 
	
	FormFront(int m) : outMode(m) {}
	
	void operator() (const IndFileLine& fline, size_t genNum) {
		generations[genNum].push_back(fline.individual);
	}
	
	void run(size_t useGenerationFirst, size_t useGenerationLast) {
		for (size_t g = useGenerationFirst; g <= useGenerationLast; ++g) {
			std::map<size_t, std::vector<GeneralOptimizationAlgorithm::Individual> >::const_iterator it;
			it = generations.find(g);
			if (it != generations.end()) {
				std::vector<size_t> dominatedCount(it->second.size(), 0);
				for (size_t i = 0; i < it->second.size(); ++i) {
					if (it->second[i].violation <= 0) {
						for (size_t j = 0; j < it->second.size(); ++j) {
							if ((it->second[i].criteria < it->second[j].criteria) && (it->second[j].violation <= 0))
								++dominatedCount[j];
						}
					}
				}
				
				for (size_t i = 0; i < it->second.size(); ++i) {
					if ((dominatedCount[i] == 0) && (it->second[i].violation <= 0))
						outputIndividual(it->second[i], g, outMode);
				}
			}
		}
	}
};


struct FindBest {
	int outMode;
	std::map<size_t, GeneralOptimizationAlgorithm::Individual> generations;
	
	FindBest(int m) : outMode(m) {}
	
	void operator() (const IndFileLine& fline, size_t genNum) {
		std::map<size_t, GeneralOptimizationAlgorithm::Individual>::iterator it = generations.find(genNum);
		if (it == generations.end())
			generations[genNum] = fline.individual;
		else if ((fline.individual.violation <= 0) && (fline.individual.criteria < it->second.criteria))
			it->second = fline.individual;
	}
	
	void run(size_t useGenerationFirst, size_t useGenerationLast) {
		for (size_t g = useGenerationFirst; g <= useGenerationLast; ++g) {
			std::map<size_t, GeneralOptimizationAlgorithm::Individual>::const_iterator it = generations.find(g);
			if (it != generations.end())
				outputIndividual(it->second, g, outMode);
		}
	}
};


struct CalcHyperVolume : public FormFront {
	CalcHyperVolume() : FormFront(out_all) {}
	std::vector<GeneralOptimizationAlgorithm::Individual>* pop;
	std::vector<double> minCrit, maxCrit;
	
	double hyperVolume(size_t gen) {
		getMinMax();
		
		std::map<size_t, std::vector<GeneralOptimizationAlgorithm::Individual> >::iterator it;
		it = generations.find(gen);
		if (it != generations.end()) 
			pop = &(it->second);
		else 
			pop = 0;
		
		if ((pop == 0) || (pop->size() == 0)) 
			return 0;
		
		HsoPopulation<GeneralOptimizationAlgorithm::Individual> hso(*pop, minCrit, maxCrit);
		return hso.getHypMeasure();
	}
	
	void getMinMax() {
		if (minCrit.size() != 0)
			return;
		
		bool uninitialized = true;
		std::map<size_t, std::vector<GeneralOptimizationAlgorithm::Individual> >::iterator it;
		for (it = generations.begin(); it != generations.end(); ++it) {
			if (uninitialized) {
				uninitialized = false;
				minCrit = maxCrit = it->second[0].criteria;
			}
			
			for (size_t i = 0; i < it->second.size(); ++i) {
				for (size_t j = 0; j < minCrit.size(); ++j) {
					if (it->second[i].violation != 0)
						continue;
					
					if (minCrit[j] > it->second[i].criteria[j]) 
						minCrit[j] = it->second[i].criteria[j];
					else if (maxCrit[j] < it->second[i].criteria[j]) 
						maxCrit[j] = it->second[i].criteria[j];
				}
			}
		}
	}
};


void runAnalysis(const std::string& individualsFile, AnalysisMode mode, const std::string params,
	int outMode, int useGenerationFirst = -1, int useGenerationLast = -1) 
{
	IndividualsFile file(individualsFile.c_str());
	
	typedef GeneralOptimizationAlgorithm::Individual Individual;
	
	std::cerr << "running analysis (mode=" << mode << ", gen=" << useGenerationFirst << "-" << 
		useGenerationLast << ")\n";
	
	{
		// analyze the file (number of generations, individuals, fields of eindivifuals)
		IndividualsFile::AnalysisResult anaRes = 
			file.analyze<Individual>(useGenerationFirst, useGenerationLast);
		if (!anaRes.bigError) {
			std::cerr << "analysis result: generations = " << anaRes.generations.front().cardinalNum 
				<< " - " << anaRes.generations.back().cardinalNum;
			if (anaRes.staticPopulationSize)
				std::cerr << ", pop. size = " << anaRes.generations.front().populationSize << "\n";
			else {
				std::cerr << ", variable pop. size:\n";
				for (size_t i = 0; i < anaRes.generations.size(); ++i) 
					std::cerr << anaRes.generations[i].populationSize << " ";
				std::cerr << "\n";
			}
		
			// did user request specific first and last generations? if yes, set
			// actual values
			if (useGenerationLast == -1) {
				useGenerationLast = anaRes.generations.back().cardinalNum;
			}
			if (useGenerationFirst == -1) {
				useGenerationFirst = anaRes.generations.back().cardinalNum;
			}
			if (useGenerationFirst > useGenerationLast)
				useGenerationLast = useGenerationFirst;
			
			if (anaRes.generations.front().cardinalNum > useGenerationFirst)
				useGenerationFirst = anaRes.generations.front().cardinalNum;
			if (anaRes.generations.back().cardinalNum < useGenerationLast)
				useGenerationLast = anaRes.generations.back().cardinalNum;
			
			IndividualsFile::FileLine<Individual> fline;
			
			// actual analysis of the file
			if (mode == an_default) {
				OutputInd outputInd(outMode);
				file.readGenerations(useGenerationFirst, useGenerationLast, fline, outputInd, anaRes.fileMode);
			} else if (mode == an_front) {
				FormFront formFront(outMode);
				file.readGenerations(useGenerationFirst, useGenerationLast, fline, formFront, anaRes.fileMode);
				formFront.run(useGenerationFirst, useGenerationLast);
			} else if (mode == an_hypervolume) {
				// hypervolume params may be specified as -params
				CalcHyperVolume hyperVol;
				if (params != "") {
					double dMin, dMax;
					char delimDummy;
					std::istringstream param(params);
					while (param && (param >> dMin >> delimDummy >> dMax)) {
						hyperVol.minCrit.push_back(dMin);
						hyperVol.maxCrit.push_back(dMax);
						param >> delimDummy;
					}
				}
				file.readGenerations(useGenerationFirst, useGenerationLast, fline, hyperVol, anaRes.fileMode);
				hyperVol.getMinMax();
				std::cout << "# hypervolume limits:\n" << "# min = [";
				for (size_t i = 0; i < hyperVol.minCrit.size(); ++i)
					std::cout << hyperVol.minCrit[i] << " ";
				std::cout << "]\n" << "# max = [";
				for (size_t i = 0; i < hyperVol.maxCrit.size(); ++i)
					std::cout << hyperVol.maxCrit[i] << " ";
				std::cout << "]\n";
				
				for (size_t gen = useGenerationFirst; gen <= (size_t)useGenerationLast; ++gen)
					std::cout << gen << " \t" << hyperVol.hyperVolume(gen) << "\n";
			} else if (mode == an_best) {
				FindBest bestSol(outMode);
				file.readGenerations(useGenerationFirst, useGenerationLast, fline, bestSol, anaRes.fileMode);
				bestSol.run(useGenerationFirst, useGenerationLast);
			}
			
		} else {
			if (anaRes.failedToOpenFile)
				std::cerr << "analysis failed (failed to open " << individualsFile << ")\n";
			else
				std::cerr << "analysis failed (file read error)\n";
		}
	}
}


int main(int argc, char** argv) {
	try {
#ifndef NO_MPI		
		Mpi::Environment mpi(argc, argv);
		Mpi::Buffer mpiBuffer(1000000);
#endif
		ArgumentsCStyle args(argc, argv);
		
		std::string settingsFname = "settings.ini";
		args.setVar(settingsFname, "-ini");
		
		std::string individuals;
		if (args.setVar(individuals, "-individuals")) {
			AnalysisMode mode = an_default;
			{
				std::string analysisMode;
				args.setVar(analysisMode, "-analysis");
				std::map<std::string, AnalysisMode> modes;
				modes["front"] = an_front;
				modes["hypervolume"] = an_hypervolume;
				modes["best"] = an_best;
				
				std::map<std::string, AnalysisMode>::const_iterator modeIt = modes.find(analysisMode);
				if (modeIt != modes.end()) {
					mode = modeIt->second;
					std::cout << "# analysis: " << analysisMode << "\n";
				}
			}
			int firstGen = -1, lastGen = -1;
			{
				std::string gen;
				if (args.setVar(gen, "-gen")) {
					std::stringstream genStream(gen);
					genStream >> firstGen;
					genStream.ignore(1);
					genStream >> lastGen;
				}
			}
			int outMode = out_all;
			{
				std::string outModeStr;
				args.setVar(outModeStr, "-out");
				std::map<std::string, OutputMode> modes;
				modes["all"] = out_all;
				modes["chromosome"] = out_chromosome;
				modes["input"] = out_chromosome;
				modes["violation"] = out_violation;
				modes["criteria"] = out_criteria;
				modes["properties"] = out_properties;
				// modes only available to 'evaluations' files
				modes["rank"] = out_rank;
				modes["evaltime"] = out_evaluation_time;
				modes["lifetime"] = out_life_time;
				
				std::map<std::string, OutputMode>::const_iterator modeIt = modes.find(outModeStr);
				if (modeIt != modes.end()) {
					std::cout << "# output: " << outModeStr << "\n";
					outMode = modeIt->second;
				}
			}
			std::string params;
			{
				args.setVar(params, "-params");
			}
			{
			    std::string temp;
				if (args.setVar(temp, "-out_no_gen")) {
                    outMode |= out_no_gen;
				}
			}
			
			runAnalysis(individuals, mode, params, outMode, firstGen, lastGen);
		} else {
			std::string temp;
			if (args.setVar(temp, "-test")) {
				UnitTests::testPDemo();
			} else {
#ifndef NO_MPI		
				runOptimizer(settingsFname, mpi);
#else
				runOptimizer(settingsFname);
#endif //NO_MPI
			}
		}
	} catch (std::exception& e) {
		std::cout << "terminal exception occured: " << e.what() << std::endl;
	} catch (const char*e) {
		std::cout << "terminal exception occured: " << e << std::endl;
	} catch(...) {
		std::cout << "unknown terminal exception caught" << std::endl;
	}
	return 0;
}
