#ifndef PARALLELNUMERICOPTIMIZER_H_INCLUDED
#define PARALLELNUMERICOPTIMIZER_H_INCLUDED


#ifndef NO_MPI
	#include "ParallelFramework.h"
#endif
#include <Random.h>
#include "Array.h"
#include "utilities.h"
#include "Individual.h"
#include "ExternalEvaluation.h"
#include "Initializer.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <sstream>
#include <cassert>
#include <fstream>
#include <stdexcept>
#include <algorithm>


#ifndef NO_MPI
template <class Chro, class Crit, class Prop>
BinaryStream& operator<< (BinaryStream& bs, const IndividualStruc<Chro, Crit, Prop>& ind) {
	return bs << ind.chromosome << ind.criteria << ind.violation << ind.properties;
}


template <class Chro, class Crit>
BinaryStream& operator<< (BinaryStream& bs, const IndividualStruc<Chro, Crit, void>& ind) {
	return bs << ind.chromosome << ind.criteria << ind.violation;
}


template <class Chro, class Crit, class Prop>
BinaryStream& operator>> (BinaryStream& bs, IndividualStruc<Chro, Crit, Prop>& ind) {
	return bs >> ind.chromosome >> ind.criteria >> ind.violation >> ind.properties;
	
}


template <class Chro, class Crit>
BinaryStream& operator>> (BinaryStream& bs, IndividualStruc<Chro, Crit, void>& ind) {
	return bs >> ind.chromosome >> ind.criteria >> ind.violation;
}
#endif


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// class FeatureFile
/// 
class FeatureFile {
	std::vector<double> workerComTimers;
	std::vector<double> workerIdleTimers;
	std::vector<std::string> workerNames;
	
public: // to be set by parent object
	int randomSeed;
	int numOfCpus;
	bool masterIsWorker;
	int masterRank;
	int queueSize;
	
public: // to be directly used by parent object
	mutable PrecisionTimerSequence evaluationSequence;
	mutable PrecisionTimer communicationTimer; // used only in parallel settings
	mutable PrecisionTimer communicationWorkerTimer; // comulative communication as seen by workers
	mutable PrecisionTimer idleTimer; 
	mutable PrecisionTimer evaluationTimer;
	mutable PrecisionTimer ioTimer;
	mutable PrecisionTimer totalTimer;
	
	static const int c_com = 1, c_eval = 2, c_idle = 3, c_io = 4, c_other = 5;
	
public:
	FeatureFile() :
		communicationTimer(evaluationSequence, c_com),
		communicationWorkerTimer(evaluationSequence, c_com),
		idleTimer(evaluationSequence, c_idle),
		evaluationTimer(evaluationSequence, c_eval),
		ioTimer(evaluationSequence, c_io),
		totalTimer(evaluationSequence, c_other) {}
	
#ifndef NO_MPI
	bool receive(int source, BinaryStream& stream) {
		std::vector<char> nodeName;
		if ((source < (int)workerComTimers.size()) && (source >= 0) && (workerComTimers[source] < 0)) {
			stream >> workerComTimers[source];
			stream >> workerIdleTimers[source];
			stream >> nodeName;
			workerNames[source].resize(nodeName.size());
			std::copy(nodeName.begin(), nodeName.end(), workerNames[source].begin());
			return true;
		} else {
			return false;
		}
	}
	
	BinaryStream sendBuffer() const {
		BinaryStream ret;
		ret << (double)communicationWorkerTimer.totalSeconds();
		ret << (double)idleTimer.totalSeconds();
		char nodeName[MPI_MAX_PROCESSOR_NAME];
		int nameLen;
		MPI_Get_processor_name(nodeName, &nameLen);
		std::vector<char> nameVec(nodeName, nodeName+nameLen);
		ret << nameVec;
		return ret;
	}
#endif

	// save function should be called by all participating CPUs
	void save(const char *fname) {
		std::ofstream file(fname);
#ifndef NO_MPI
		try {
            // take the default communicator
            Mpi::Communicator comm;
			Gatherer gatherer(comm);
			if (comm.getRank() == 0) {
				file << "# CPUs = " << numOfCpus << "\n";
				file << "queue size (per CPU) = " << queueSize << "\n";
				file << "random seed = " << randomSeed << "\n";
				
				workerComTimers.resize(comm.getSize(), -1.0);
				workerIdleTimers.resize(comm.getSize(), -1.0);
				workerNames.resize(comm.getSize());
				
				// gather data from workers (for now this is timer data only)
				while (!gatherer.gatherIn(*this, *this, 0)) {};
				// gather data from master (itself)
				{
					BinaryStream temp = sendBuffer();
					receive(0, temp);
				}
				if (!masterIsWorker)
					file << "master: " << masterRank << " ... " << workerNames[masterRank] << "\n";
				file << "workers:\n";
				for (size_t i = 0; i < workerNames.size(); ++i) {
					if ((int)i == masterRank) {
						if (masterIsWorker)
							file << "   " << i << " ... " << workerNames[i] << " (also master)\n";
					} else
						file << "   " << i << " ... " << workerNames[i] << "\n";
				}
				
				double workerComTime = 0.0f;
				double workerIdleTime = 0.0f;
				for (size_t i = 1; i < workerComTimers.size(); ++i) {
					workerComTime += workerComTimers[i];
					workerIdleTime += workerIdleTimers[i];
				}
				double comTime = communicationTimer.totalSeconds();
				double evalTime = evaluationTimer.totalSeconds();
				double idleTime = idleTimer.totalSeconds();
				double ioTime = ioTimer.totalSeconds();
				double totalTime = totalTimer.totalSeconds();
				double otherTime = totalTime - ioTime - evalTime - comTime - idleTime;
				file << 
					"timer analysis:\n" << 
					"   communication: " << comTime << "\n" <<
					"   worker comm.:  " << workerComTime << "\n" << 
					"   idle (waiting):" << idleTime << "\n" << 
					"   worker idle:   " << workerIdleTime << "\n" << 
					"   evaluation:    " << evalTime << "\n" << 
					"   in/out:        " << ioTime << "\n" <<
					"   other:         " << otherTime << "\n" << 
					"  ---------------------------\n" <<
					"   total:         " << totalTime << "\n";
			} else {
				gatherer.gatherIn(*this, *this, 0);
			}
		} catch (const char* e) {
			std::cout << "Error while gathering worker info:\n";
			std::cout << "   " << e << "\n";
		} catch (...) {
			std::cout << "Unknown error while gathering worker info:\n";
		}
#else
		file << "# CPUs = 1, MPI disabled\n";
		file << "random seed = " << randomSeed << "\n";

		double comTime = communicationTimer.totalSeconds(); // must be 0
		double evalTime = evaluationTimer.totalSeconds() - comTime;
		double ioTime = ioTimer.totalSeconds();
		double totalTime = totalTimer.totalSeconds();
		double otherTime = totalTime - ioTime - evalTime - comTime;
		file << "timer analysis:\n" << 
			"   evaluation:    " << evalTime << "\n" << 
			"   in/out:        " << ioTime << "\n" <<
			"   other:         " << otherTime << "\n" << 
			"  ---------------------------\n" <<
			"   total:         " << totalTime << "\n";
#endif
        file << std::flush;
	}
};


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// template class ParallelNumericOptimizer
/// 
/// Float:
///  - a simple numeric class in which trial solutions are encoded and algebraic 
/// 	operations (a reasonable choice for this class is the built-in double)
/// ChromosomeSize:
///  - number of Float genes in a chromosome (number of parameters representing a 
/// 	solution)
/// Algorithm:
///  - class that encapsulates the inner workings of an algorithm
///  - 
/// Evaluation:
///  - provides "operator(const Chromosome&)" through which evaluation of a trial 
/// 	solution is executed
///  - provides "class Value", that defines the return value of "operator()"
///  - provides function "normalize(Chromosome& solution) const" for normalization 
///		of individuals (keeping their genes within predefined bounds)
/// 
template<class Algorithm, class Chromosome, class Evaluation>
class ParallelNumericOptimizer {
protected:
	typedef typename Evaluation::Value FunctionValue;
	typedef typename Evaluation::Properties Properties;
	
public: // to be used by different implementations of NumericOptimizer
	typedef IndividualStruc<Chromosome, FunctionValue, Properties> Individual;
	
protected:
	std::vector<Individual> pop, candidates;
	Evaluation eval;
	size_t currentEvaluation;
	size_t totalGenerated;
	size_t numTruncations;
	bool (ParallelNumericOptimizer::*terminationCriteriaMet)() const;
	
public: // to be set directly from the parent object
	size_t maxNumOfEvaluations;
	size_t evaluationsPerTruncation;
	size_t populationSize;
	// here provided random seed is used by numeric optimizer at the population 
	// initialization (first thing to possibly use it)
	// there should be no changes of the random generator used in optimizer (currently
	// this is CRand) between the calls to initialize population and evolution
	int randomSeed;
	int minQueueLength, maxQueueLength;
	
	FeatureFile featureFile;
	bool debugOutput;
	
public:
	ParallelNumericOptimizer();
	
	template<class Initializer>
	void initPopulation(size_t popSize, const Initializer& input);
#ifndef NO_MPI
	void evolution(const Mpi::Communicator& com = MPI_COMM_WORLD);
#else
	void evolution();
#endif //NO_MPI
	
	const std::vector<Individual>& population() const {return pop;};
	Evaluation& evaluation() {return eval;}
	
protected:
	// output
	void outputIndividuals() const;
	
	// termination criteria
	bool limitedNumOfEvaluations() const;
	
	// members for parallel framework:
	std::map<int, int> indexIdMap;
	
public:
	Chromosome generateJob(int id);
	int moreJobs(int preferredNum);
	// when return value is true, old jobs can be invalidated
	bool newResult(Individual& out, int id, int rank, double evolTime, double lifeTime);
	
protected:
	struct EvaluationJob {
		ParallelNumericOptimizer& parent;
		
		typedef typename Evaluation::Input Input;
		typedef typename Evaluation::Value Output;
		typedef typename Evaluation::Properties Properties;
		typedef IndividualStruc<Input, Output, Properties> Individual;
		
		EvaluationJob(ParallelNumericOptimizer& p) : parent(p) {}
		
		inline Individual operator()(Input& el) {
			try {
				ScopeTimer t(parent.featureFile.evaluationTimer);
				Individual out;
				parent.eval.normalize(el);
				if (Individual::hasProperties)
					parent.eval(el, out.criteria, out.violation, out.getProperties());
				else
					parent.eval(el, out.criteria, out.violation);
				out.chromosome = el;
				return out;
			} catch (std::exception& e) {
				std::cout << "\n exception during evaluation:\n " << e.what() << std::endl;
				throw;
			} catch (...) {
				std::cout << "\n unknown exception during evaluation " << std::endl;
				throw;
			}
		}
	};
};


// #################################################################################
// implementation below
// #################################################################################


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// template class ParallelNumericOptimizer
template<class Algorithm, class Chromosome, class Evaluation>
ParallelNumericOptimizer<Algorithm, Chromosome, Evaluation>::ParallelNumericOptimizer() {
	terminationCriteriaMet = &ParallelNumericOptimizer::limitedNumOfEvaluations;

	currentEvaluation = 0;
	totalGenerated = 0;
	numTruncations = 0;
	
	maxNumOfEvaluations = 1000;
	evaluationsPerTruncation = populationSize = 10;
	minQueueLength = maxQueueLength = 1;
}


template<class Algorithm, class Chromosome, class Evaluation> template<class Initializer>
void ParallelNumericOptimizer<Algorithm, Chromosome, Evaluation>::initPopulation(size_t popSize, const Initializer& input) {
	//ScopeTimer t(featureFile.totalTimer);
	using namespace Random;
	
	// set random seed
	if (randomSeed == 0) {
		randomSeed = CRand::randomizeSeed();
	} else {
		CRand::setSeed(randomSeed);
	} 
	featureFile.randomSeed = randomSeed;
	
	evaluationsPerTruncation = populationSize = popSize;
	pop.resize(populationSize);
	for (size_t i = 0; i < pop.size(); ++i)
		input.generateSolution(pop[i]);
}


template<class Algorithm, class Chromosome, class Evaluation>
#ifndef NO_MPI
void ParallelNumericOptimizer<Algorithm, Chromosome, Evaluation>::evolution(const Mpi::Communicator& com) {
#else
void ParallelNumericOptimizer<Algorithm, Chromosome, Evaluation>::evolution() {
#endif //NO_MPI
	{
		ScopeTimer t(featureFile.totalTimer);
		EvaluationJob job(*this);
		
#ifndef NO_MPI
		PFramework<Chromosome, Individual> fw(com);
		fw.debugOutput = debugOutput;
		featureFile.numOfCpus = com.getSize();
		featureFile.masterIsWorker = fw.masterIsWorker;
		featureFile.masterRank = fw.masterRank;
		featureFile.queueSize = minQueueLength;
		fw.setQueueSizes(minQueueLength, maxQueueLength);
		
		// algorithm main loop
		while (fw.execute(job, *this, (com.getRank() == 0 ? &featureFile.communicationTimer : 
			&featureFile.communicationWorkerTimer), &featureFile.idleTimer))
		{
			#ifdef WIN32
			// used when testing on single CPU windows machine - Sleep(0) yields current 
			// process CPU time
			Sleep(0);
			#endif
		}
#else
		featureFile.numOfCpus = 1;
		featureFile.masterIsWorker = true;
		featureFile.masterRank = 0;
		featureFile.queueSize = minQueueLength;
		while (moreJobs(1) > 0) {
			static int lastId = 0;
			lastId++;
			// generate & execute job
			Chromosome c = generateJob(lastId);
			PrecisionTimer t;
			t.start();
			Individual ind = job(c);
			// register the result
			t.pause();
			newResult(ind, lastId, 0, t.totalSeconds(), 0);
		}
#endif //~NO_MPI
	}
	
#ifndef NO_MPI
    // this code is not requ ired for funcionallity of AMS-DEMO and is not timed (it could go under iotimer). 
    {
        std::ofstream seqFile;
        std::ostringstream seqName;
        seqName << "timeSeq" << com.getRank() << ".txt";
        seqFile.open(seqName.str().c_str());
        for (size_t i = 0; i < featureFile.evaluationSequence.size(); ++i) {
            seqFile << std::fixed << std::setprecision(6) << featureFile.evaluationSequence.at(i).first;
            seqFile << " " << featureFile.evaluationSequence.at(i).second << "\n";
        }
    }
#endif //NO_MPI
	
	// save feature file (including gathering timer infos from workers)
	featureFile.save("featureFile.txt");
}


template<class Algorithm, class Chromosome, class Evaluation>
void ParallelNumericOptimizer<Algorithm, Chromosome, Evaluation>::outputIndividuals() const {
	ScopeTimer t(featureFile.ioTimer);
	
	IndividualsFile file("individuals.txt");
	file.outputGeneration(pop, numTruncations);
}


template<class Algorithm, class Chromosome, class Evaluation>
bool ParallelNumericOptimizer<Algorithm, Chromosome, Evaluation>::limitedNumOfEvaluations() const {
	return (currentEvaluation >= maxNumOfEvaluations);
}


template<class Algorithm, class Chromosome, class Evaluation>
Chromosome ParallelNumericOptimizer<Algorithm, Chromosome, Evaluation>::generateJob(int id) {
//	std::cerr << "generating candidate " << id << " (" << totalGenerated << ")\n";
	if (totalGenerated < populationSize) {
		// first population
		++totalGenerated;
		indexIdMap[id] = totalGenerated-1;
//		std::cerr << "generated " << indexIdMap[id] << "\n";
//		std::cerr << "generated candidate from 1st pop: " << pop[totalGenerated-1].chromosome << std::endl;
		return pop[totalGenerated-1].chromosome;
	} else {
		++totalGenerated;
		Chromosome c = static_cast<Algorithm*>(this)->generateCandidate();
//		std::cerr << "generated candidate " << c << std::endl;

//		static std::ofstream prva("prva.txt", std::ios::app);
//		for (size_t i = 0; i < c.size(); ++i) {
//			if (!std::isfinite(c[i])) {
//				prva << "generated " << c << std::endl;
//				break;
//			}
//		}
		return c;
	}
}


template<class Algorithm, class Chromosome, class Evaluation>
int ParallelNumericOptimizer<Algorithm, Chromosome, Evaluation>::moreJobs(int preferredNum) {
	if (currentEvaluation >= maxNumOfEvaluations)
		return -1;
	return std::min(preferredNum, (int)maxNumOfEvaluations-(int)totalGenerated);
}
	

template<class Algorithm, class Chromosome, class Evaluation>
bool ParallelNumericOptimizer<Algorithm, Chromosome, Evaluation>::newResult(Individual& out, int id, int rank, double evolTime, double lifeTime) {
	// !!! evaluirana populacija se ne shranjuje pravilno nazaj v populacijo
//	std::cerr << "receiving result " << id << " (" << currentEvaluation << ")\n";
//	std::cerr << "receiving result " << out.chromosome << " | " << out.criteria << ")\n";
	bool nextGeneration = false;
	{
        ScopeTimer t(featureFile.ioTimer);
		IndividualsFile file("evaluations.txt");
		file.outputEvaluation(out, currentEvaluation, rank, evolTime, lifeTime);
	}
	
	currentEvaluation++;
	std::map<int, int>::iterator indexIt = indexIdMap.find(id);
	if (indexIt != indexIdMap.end()) {
		// first generation -> just save the normalized & evaluated individual
		// back into the population

		pop[indexIt->second] = out;
		indexIdMap.erase(indexIt);
	} else {
		// parent index should be stored alongside each trial individual (in a similar way
		// as is already implemented for the first generation);
		// after each truncation (and shuffeling of individuals), the stored parent indices
		// should be rechecked
		size_t parentIndex = (currentEvaluation-1+indexIdMap.size()) % populationSize;
		std::map<int, int>::iterator parentIt = indexIdMap.find(parentIndex);
		
		if (pop[parentIndex].evaluated()) {
			static_cast<Algorithm*>(this)->performSelection(parentIndex, out);
		} else {
			// problematic handling of unevaluated first generation individuals
			// replace parent with evaluated individual and hope parent will
			// be evaluated soon and trialed against this individual
			// (otherwise population diversity will decrease)
			
			for (std::map<int, int>::iterator it = indexIdMap.begin(); it != indexIdMap.end(); ++it) {
				if (it->first == (int)parentIndex) {
					indexIdMap.erase(it);
					break;
				}
			}
			pop[parentIndex] = out;
		}
	}
	
	if (currentEvaluation % evaluationsPerTruncation == 0) {
		if (currentEvaluation > evaluationsPerTruncation) {
			nextGeneration = true;
			++numTruncations;
			static_cast<Algorithm*>(this)->performTruncate();
		}
		outputIndividuals();
	}
	
	return nextGeneration;
}


#endif // PARALLELNUMERICOPTIMIZER_H_INCLUDED
