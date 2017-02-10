#include "sim.h"


#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <cctype>
#include <ctime>


#include <Ini.h>
#include <Arguments.h>
#include "vectorMath.h"
#include "nonlinearFit.h"
#include "columnFile.h"
#include <AMS-DEMO/GeneralOptimizationAlgorithm.cpp>
#include <AMS-DEMO/utilities.cpp>


#define DBVAL(A) std::cout << #A << " = " << (A) << " ";
#define DBVALn(A) DBVAL(A); std::cout << "\n";


using namespace std;


/*
typedef std::vector<double> VecDouble;
typedef std::vector<TFunction> SimResult;
*/

#include "sim.h"
#include "simlib/sim_lib.h"
#include "nonlinearFit.h"


void test() {
	testWohlInterpolation();
}


/**
    Here, all the arguments and settings read from *.ini are specified
**/
struct ArgumentSystem {
	struct mapElement {
		std::string description;
		void (ArgumentSystem::*function)();
	};

	// settings read through arguments
	struct Settings {
		enum Mode {
			mode_help = 1,
			mode_test = 2,
			mode_single_sim = 4,
			mode_optimization = 8,
			mode_unknown = -1
		};

		std::string iniFname;
		std::string freeFname;
		int mode;
		std::vector<double> simParams;

		OutputSettings outSettings;

	public:
		Settings() :
			iniFname("settings.ini"),
			freeFname(""),
			mode(8)
		{
		}
	} settings;

	DashArguments args;
	std::map<std::string, mapElement> argsMap;
	typedef std::map<std::string, mapElement>::iterator mapIter;
	typedef std::map<std::string, mapElement>::iterator constMapIter;
	std::vector<std::string> activeParams;

public:
    /// setup DashArguments
	ArgumentSystem(int argc, char** argv) : args(argc, argv) {
		add("-test",
			"run test function (changable function used for debugging)",
			&ArgumentSystem::setTest);
		add("-?",
			"show this help screen",
			&ArgumentSystem::setHelp);
		add("-sim",
			"just run single a simulation with parameters provided after -sim\n",
			&ArgumentSystem::setSim);
		add("-out",
			"specify outputs of the program; possible values include result, "
			"layer_aps, cell_aps <num> [<num>]*, excitation_sequence"
			"multiple outputs are allowed, e.g. -out result -out cell_aps 1",
			&ArgumentSystem::setOutputs);
	}

	/// output help to std.out
	void showHelp() {
		std::cout << "Argument list:\n";
		for (mapIter it = argsMap.begin(); it != argsMap.end(); ++it) {
			std::cout << "   " << it->first << " \t" << it->second.description << "\n";
		}
	}

	/// parse arguments that are not a part of the DashArguments
	void parse() {
		std::cerr << "***** parsing program arguments ****************************\n";
		// first free argument is considered a filename (not yet clear what to do with this file)
		if (args.freeArguments().size() > 0)
		settings.freeFname = args.freeArguments()[0];
		//DBVALn(settings.freeFname);

		for (constMapIter it = argsMap.begin(); it != argsMap.end(); ++it) {
			if (args.isSet(it->first)) {
				activeParams.clear();
				args.setVars(activeParams, it->first);
				(this->*it->second.function)();
			}
		}

		if (((settings.mode & Settings::mode_help) > 0) || (settings.mode == 0)){
			// help is of highest priority, only show help, and quit
			std::cerr << " displaying help screen\n";
			showHelp();
			std::cerr << "\n";
			settings.mode ^= Settings::mode_help;
		}
		if ((settings.mode & Settings::mode_test) > 0) {
			// second priority, run test function and quit
			std::cerr << " running test function\n";
			test();
			settings.mode ^= Settings::mode_test;
		} else {
			if ((settings.mode & Settings::mode_single_sim) > 0) {
				settings.mode &= ~Settings::mode_optimization;

				if (args.isSet("-sim")) {
					std::istringstream simArgStr(args["-sim"]);
					settings.simParams.clear();
					for (double temp; simArgStr >> temp; ) {
						settings.simParams.push_back(temp);
						char tempch = simArgStr.peek();
						if ((tempch == ',') || (tempch == ';'))
							simArgStr.ignore(1);
					}
					std::cout << " parameters for single simulator run: " << settings.simParams << "\n";
                    if (settings.simParams.size() < 1) {
                        throw std::runtime_error(" Error: not enough parameters to run simulation: "+args["-sim"]);
                    }
				} else
					throw std::runtime_error(" Error: single sim mode but unknown settings");
			} else {
			}
		}

		std::string individuals ("individuals.txt");
		if (args.setVar(individuals, "-individuals")) {
		}
		std::cerr << "\n";
	}

private:
    /// specify a non-dash argument
	void add(const std::string& name, const std::string& description, void (ArgumentSystem::*function)()) {
		mapElement temp;
		temp.description = description;
		temp.function = function;
		argsMap.insert(std::make_pair(name, temp));
	}

	void setTest() 		{ settings.mode |= Settings::mode_test; }
	void setSim() 		{ settings.mode |= Settings::mode_single_sim; }
	void setHelp()		{ settings.mode |= Settings::mode_help; }
	void setOutputs() {
		for (std::vector<std::string>::const_iterator it = activeParams.begin();
			it != activeParams.end();
			++it)
		{
			std::istringstream stream(*it);
			std::string outputName;
			stream >> outputName;
			//std::cout << outputName << "\n";
			if (outputName == "cell_aps") {
				while (stream) {
					stream.ignore(1);
					size_t num;
					stream >> num;
					if (stream)
						settings.outSettings.outputCellAps.push_back(num);
				}
			} else if (outputName == "layer_aps") {
				settings.outSettings.layerAps = true;
//			} else if (outputName == "excitation_sequence") {
//				settings.outSettings.excitationSequence = true;
			} else if (outputName == "result") {
				settings.outSettings.result = true;
			} else {
				std::cerr << "skipping an unrecognized output opition: " << outputName << "\n";
			}
		}
		std::cout << "\n";
	}
};


/**
    First of the two main functions. If optimization is requested, this function is run.
**/
void runOptimization(int argc, char** argv) {
	int randSeed = 0;
#ifndef NO_MPI
	Mpi::Environment mpi(argc, argv);
	Mpi::Buffer mpiBuffer(100000);
	{
		char nodeName[MPI_MAX_PROCESSOR_NAME];
		int nameLen, nameSum = 0;
		MPI_Get_processor_name(nodeName, &nameLen);
		for (int i=0; i < nameLen; ++i)
			nameSum += ((int)nodeName[i]) << (i % 24);
		randSeed = time(NULL) + nameSum;
	}
#endif
	Random::CRand::randomizeSeed();
	GeneralOptimizationAlgorithm alg;
	{
		alg.settings.evaluation.propertiesVectorLength = 0;
		std::shared_ptr<OptimizationFunction> optFunct(
			new OptimizationFunction(alg.settings.evaluation.propertiesVectorLength));
		optFunct->getGeneParams(
				alg.settings.evaluation.chromosomeVectorLength,
				alg.settings.evaluation.criteriaVectorLength,
				alg.settings.initialPopulation.geneMin,
				alg.settings.initialPopulation.geneMax);

		alg.settings.algorithmParams.randomSeed = randSeed;
		optFunct->getEvolutionaryParams(
			alg.settings.algorithmParams.maxNumOfEvaluations,
			alg.settings.algorithmParams.populationSize,
			alg.settings.algorithmParams.minQueueLength);
		if (alg.settings.algorithmParams.populationSize == 0)
			alg.settings.algorithmParams.populationSize =
				int(floor(alg.settings.evaluation.chromosomeVectorLength
				* sqrt(double(alg.settings.evaluation.criteriaVectorLength)) * 2.0));
		alg.settings.algorithmParams.maxQueueLength = alg.settings.algorithmParams.minQueueLength;
		alg.settings.algorithmParams.maxNumOfEvaluations *=
			alg.settings.algorithmParams.populationSize;
		alg.settings.algorithmParams.crossoverProbability = 0.5;
		alg.settings.algorithmParams.scalingFactors.push_back(0.5);
		alg.settings.algorithmParams.deSchema = "rand/1/bin";

		alg.settings.evaluation.internalFunction = optFunct;
        alg.settings.evaluation.internalFunctionName = "EkgSim";
	}

	if (alg.updateSettings()) {
#ifdef NO_MPI
		std::cerr << "##### Running EkgSim v3 based optimization, "
			<< (alg.settings.evaluation.criteriaVectorLength == 1 ? "DE ##" : "DEMO " )
			<< "#################\n" << std::flush;

		// run calls:
		//	- population initialization
		//	- evolution
		alg.run();
#else
		alg.setCommunicator(mpi.getCommunicator());

		if (mpi.getCommunicator().getRank() == 0) {
			std::cerr << "##### Running EkgSim v3 based optimization, "
				<< (alg.settings.evaluation.criteriaVectorLength == 1 ? "DE" : "DEMO" )
				<< " with MPI on " << mpi.getCommunicator().getSize() << " CPUs #####\n"
				<< std::flush;
			alg.run();
		} else {
			alg.runWorker();
		}
#endif
	} else {
		std::cout << "error in settings, could not setup evolutionary algorithm\n";
	}
}


/**
    Second of the two main functions. If simulation is requested and parameters are given, this function is run.
**/
void runSingleSim(const ArgumentSystem::Settings& settings) {
	std::cerr << "##### Running a single simulation experiment ####################\n";

	// declare optimization function from sim.h
	OptimizationFunction func(0);
	func.setup(settings.outSettings);

	std::vector<double> result, params;

	double violation;
	PrecisionTimer pt;
	{
		ScopeTimer t(pt);
		func(settings.simParams, result, violation, params);
	}
	std::cout << " simulation done in " << pt.totalSeconds() << " seconds\n";
	std::cout << " criteria = " << result << ", violation = " << violation << "\n";
}


int main(int argc, char** argv) {
	try {
		ArgumentSystem arg(argc, argv);
		arg.parse();

		switch (arg.settings.mode) {
			case ArgumentSystem::Settings::mode_single_sim:
				runSingleSim(arg.settings);
				break;
			case ArgumentSystem::Settings::mode_optimization:
				runOptimization(argc, argv);
				break;
			default:
				std::cout << "missing program arguments, stopping execution\n";
				break;
		}
        
        std::cout << "All done\n";
	} catch (const char* e) {
		std::cout << "exception caught: " << e << "\n";
    } catch (std::runtime_error& e) {
		std::cout << "runtime error caught: " << e.what() << "\n";
	} catch (std::exception& e) {
		std::cout << "std exception caught: " << e.what() << "\n";
	} catch (...) {
		std::cout << "unknown exception caught, execution halted\n";
	}
}
