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

#ifndef NUMERICOPTIMIZER_H_INCLUDED
#define NUMERICOPTIMIZER_H_INCLUDED


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
#include "Random.h"
#include "Array.h"
#include "utilities.h"
#include "Individual.h"
#include "ExternalEvaluation.h"
#include "Initializer.h"
#ifndef NO_MPI
	#include "jobDistributer.h"
#endif //~NO_MPI


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// template class NumericOptimizer
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
class NumericOptimizer {
protected:
	typedef typename Evaluation::Value FunctionValue;
	
public:
	typedef IndividualStruc<Chromosome, FunctionValue> Individual;
	
protected:
	std::vector<Individual> pop, children;
	Evaluation eval;
	size_t currentGeneration;
	bool (NumericOptimizer::*terminationCriteriaMet)() const;
	
public:
	size_t maxNumOfGenerations;
	
	// comulative timers for different parts of algorithm
	mutable PrecisionTimerSequence evaluationSequence;
	
	mutable PrecisionTimer communicationTimer; // used only in parallel settings
	mutable PrecisionTimer communicationWorkerTimer; // comulative communication as seen by workers
	mutable PrecisionTimer evaluationTimer;
	mutable PrecisionTimer ioTimer;
	mutable PrecisionTimer totalTimer;
	
	static const int c_com = 1, c_eval = 2, c_io = 4, c_other = 5;
	
	
public:
	NumericOptimizer();
	
	template<class Initializer>
	void initPopulation(size_t populationSize, const Initializer& input);
	void evolution();
#ifndef NO_MPI
	void parallelEvolution(const Mpi::Communicator& com, bool workingMaster = false);
#endif
	
	const std::vector<Individual>& population() const {return pop;};
	Evaluation& evaluation() {return eval;}
	
protected:
	void evaluatePopulation(std::vector<Individual>& pop);
#ifndef NO_MPI
	void parallelEvaluatePopulation(std::vector<Individual>& pop, const Mpi::Communicator& com, bool workingMaster);
#endif	
	void outputIndividuals() const;
	
	// termination criteria
	bool limitedNumOfGenerations() const;
};


// #################################################################################
// implementation below
// #################################################################################


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// template class NumericOptimizer
template<class Algorithm, class Chromosome, class Evaluation>
NumericOptimizer<Algorithm, Chromosome, Evaluation>::NumericOptimizer() : 
	communicationTimer(evaluationSequence, c_com),
	communicationWorkerTimer(evaluationSequence, c_com),
	evaluationTimer(evaluationSequence, c_eval),
	ioTimer(evaluationSequence, c_io),
	totalTimer(evaluationSequence, c_other)
{
	terminationCriteriaMet = &NumericOptimizer::limitedNumOfGenerations;
	maxNumOfGenerations = 100;
}


template<class Algorithm, class Chromosome, class Evaluation> template<class Initializer>
void NumericOptimizer<Algorithm, Chromosome, Evaluation>::initPopulation(size_t populationSize, const Initializer& input) {
	ScopeTimer t(totalTimer);
	
	pop.resize(populationSize);
	children.resize(populationSize);
	for (size_t i = 0; i < pop.size(); ++i)
		input.generateSolution(pop[i]);
	children = pop;
}


template<class Algorithm, class Chromosome, class Evaluation>
void NumericOptimizer<Algorithm, Chromosome, Evaluation>::evolution() {
	{
		ScopeTimer t(totalTimer);
	
		// step through generations
		evaluatePopulation(pop);
		for(this->currentGeneration = 0; !(this->*this->terminationCriteriaMet)(); ++this->currentGeneration) { 
			outputIndividuals();
			// generate new generation and put it into "children"
			static_cast<Algorithm*>(this)->nextGeneration();
			evaluatePopulation(children);
			// perform selection and fill population "pop" with most fit individuals from 
			// old population "pop" and new population "children"
			static_cast<Algorithm*>(this)->generationalSelection();
		}
		outputIndividuals();
	}

	{	
		float evalTime = evaluationTimer.totalSeconds();
		float ioTime = ioTimer.totalSeconds();
		float totalTime = totalTimer.totalSeconds();
		float otherTime = totalTime - ioTime - evalTime;
		std::cerr << "timer analysis:\n" << 
			"   evaluation:    " << evalTime << "\n" << 
			"   in/out:        " << ioTime << "\n" <<
			"   other:         " << otherTime << "\n" << 
			"  ---------------------------\n" <<
			"   total:         " << totalTime << "\n";
	}
}


template<class Algorithm, class Chromosome, class Evaluation>
void NumericOptimizer<Algorithm, Chromosome, Evaluation>::evaluatePopulation(std::vector<Individual>& p) {
	ScopeTimer t(evaluationTimer);
	for (size_t i = 0; i < p.size(); ++i) {
		if (p[i].violation == -1)
			p[i].violation = eval(p[i].chromosome, p[i].criteria);
	}
}


#ifndef NO_MPI
namespace Hidden {
	template<class Evaluation>
	struct Job {
		Evaluation& eval;
		
		typedef typename Evaluation::Input Input;
		typedef typename Evaluation::Value Output;
		
		Job(Evaluation& e) : eval(e) {}
		
		inline double operator()(const Input& el, Output& out) {
			try {
				//std::cout << el << std::endl;
				return eval(el, out);
			} catch (std::exception& e) {
				std::cout << "\n exception during evaluation:\n " << e.what() << std::endl;
				throw;
			} catch (...) {
				std::cout << "\n unknown exception during evaluation " << std::endl;
				throw;
			}
		}
	};


	template<class Algorithm, class Chromosome, class Evaluation>
	struct WorkerFeedbackJob {
		typedef Array<float, 1> Input;
		typedef Array<float, 1> Output;
		
		NumericOptimizer<Algorithm, Chromosome, Evaluation>& parent;
		
		WorkerFeedbackJob(NumericOptimizer<Algorithm, Chromosome, Evaluation>& p) : parent(p) {}
		
		inline double operator()(const Input&, Output& out) {
			out[0] = parent.communicationWorkerTimer.totalSeconds();
			return 0.0;
		}
	};

}


template<class Algorithm, class Chromosome, class Evaluation>
void NumericOptimizer<Algorithm, Chromosome, Evaluation>::parallelEvolution(const Mpi::Communicator& com, bool workingMaster) {
	using namespace Hidden;
	{
		ScopeTimer t(totalTimer);
		
		// workingMaster should be true if master process is doing work too; if false, master 
		// process only delegates work to other processes and then waits for results
		
		// step through generations
		parallelEvaluatePopulation(pop, com, workingMaster);
		if (com.getRank() == 0) {
			for(this->currentGeneration = 0; !(this->*this->terminationCriteriaMet)(); ++this->currentGeneration) { 
				outputIndividuals();
				// generate new generation and put it into "children"
				static_cast<Algorithm*>(this)->nextGeneration();
				
				parallelEvaluatePopulation(children, com, workingMaster);
				// perform selection and fill population "pop" with most fit individuals from 
				// old population "pop" and new population "children"
				static_cast<Algorithm*>(this)->generationalSelection();
			}
			outputIndividuals();
		} else {
			for(this->currentGeneration = 0; !(this->*this->terminationCriteriaMet)(); ++this->currentGeneration)
				parallelEvaluatePopulation(children, com, workingMaster);
		}
	}
	
	// gather timer infos from workers and output all timer infos
	{	
		JobDistributer<Array<float, 1>, Array<float, 1> > distributer(com);
		WorkerFeedbackJob<Algorithm, Chromosome, Evaluation> job(*this);
		
		if (distributer.comm().getRank() == 0) {
			std::vector<IndividualStruc<Array<float, 1>, Array<float, 1> > > times;
			times.resize(distributer.comm().getSize());
			for (size_t i = 1; i < times.size(); ++i)
				times[i].violation = -1;
			
			distributer.executeBoss(times, job);
			
			float workerComTime =0.0f;
			for (size_t i = 1; i < times.size(); ++i) 
				workerComTime += times[i].criteria[0];
			float comTime = communicationTimer.totalSeconds();
			float evalTime = evaluationTimer.totalSeconds() - comTime;
			float ioTime = ioTimer.totalSeconds();
			float totalTime = totalTimer.totalSeconds();
			float otherTime = totalTime - ioTime - evalTime - comTime;
			std::cerr << "timer analysis:\n" << 
				"   communication: " << comTime << "\n" <<
				"   worker comm.:  " << workerComTime << "\n" << 
				"   evaluation:    " << evalTime << "\n" << 
				"   in/out:        " << ioTime << "\n" <<
				"   other:         " << otherTime << "\n" << 
				"  ---------------------------\n" <<
				"   total:         " << totalTime << "\n";
		} else {
			distributer.executeWorker(job);
			/*
			float workerComTime = communicationWorkerTimer.totalSeconds();
			float evalTime = evaluationTimer.totalSeconds() - workerComTime;
			
			std::cerr << "worker " << distributer.comm().getRank() << " analysis:\n" <<
				"   worker comm.:  " << workerComTime << "\n" << 
				"   evaluation:    " << evalTime << "\n";
				*/
		}
	}
	{
		std::ofstream seqFile;
		std::ostringstream seqName;
		seqName << "timeSeq" << com.getRank() << ".txt";
		seqFile.open(seqName.str().c_str());
		for (size_t i = 0; i < evaluationSequence.size(); ++i) {
			seqFile << std::fixed << std::setprecision(6) << evaluationSequence.at(i).first;
			seqFile << " " << evaluationSequence.at(i).second << "\n";
		}
	}
}


template<class Algorithm, class Chromosome, class Evaluation>
void NumericOptimizer<Algorithm, Chromosome, Evaluation>::parallelEvaluatePopulation(std::vector<Individual>& p, const Mpi::Communicator& com, bool workingMaster) {
	using namespace Hidden;
	
	ScopeTimer t(evaluationTimer);
	
	JobDistributer<Chromosome, FunctionValue> distributer(com);
	Job<Evaluation> job(eval);
	
	if (distributer.comm().getRank() == 0) {
		if (workingMaster)
			distributer.executeBossAndWorker(p, job, &communicationTimer);
		else
			distributer.executeBoss(p, job, &communicationTimer);
	} else {
		distributer.executeWorker(job, &communicationWorkerTimer);
	}
}
#endif //NO_MPI


template<class Algorithm, class Chromosome, class Evaluation>
void NumericOptimizer<Algorithm, Chromosome, Evaluation>::outputIndividuals() const {
	ScopeTimer t(ioTimer);
	IndividualsFile file("individuals.txt");
	file.outputGeneration(pop, currentGeneration);
}


template<class Algorithm, class Chromosome, class Evaluation>
bool NumericOptimizer<Algorithm, Chromosome, Evaluation>::limitedNumOfGenerations() const {
	return (currentGeneration >= maxNumOfGenerations);
}


#endif // NUMERICOPTIMIZER_H_INCLUDED
