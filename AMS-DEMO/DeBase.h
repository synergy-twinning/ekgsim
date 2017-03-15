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

#ifndef DEBASE_H_INCLUDED
#define DEBASE_H_INCLUDED


#include <iostream>
#include <vector>
#include <map>
#include "Array.h"
#include <Random.h>
#include "Individual.h"


template <class Chromosome>
struct Element {
	typedef float Type;
};

template <class Chromosome, class FunctVal>
struct Element<IndividualStruc<Chromosome, FunctVal> > {
	typedef typename Element<Chromosome>::Type Type;
};

template <class Float, size_t N>
struct Element<Array<Float, N> > {
	typedef Float Type;
};


/// ##################################################################################
/// template class DeBase
/// 
/// definition of DE parameters:
///  - scaling factors (as vector)
///  - crossover probability
///
/// definition of DE help structures:
///  - indices (used to randomize population order)
///
/// implementation of trial solution generation:
///  - selection of base individual (parent)
///  - generation of candidate individual (mutation)
///  - recombination of parent and trial individuals
///
template<class Individual, class Rand = Random::CRand>
class DeBase {
protected:
	typedef void (DeBase::*SelectBase)(size_t& baseIndex) const;
	typedef void (DeBase::*GenerateCandidate)(const std::vector<Individual>& pop, size_t parentIndex, Individual& dest) const;
	typedef void (DeBase::*RecombineCandidate)(const Individual& parent, Individual& dest) const;

	SelectBase selectBase;
	GenerateCandidate generateCandidate;
	RecombineCandidate recombineCandidate;
	
	std::vector<size_t> indices;
	std::string schema;
	
	typedef typename Element<Individual>::Type Float;
	
public:
	// algorithm parameters
	std::vector<Float> scalingFactor;
	double crossoverProbability;
	
public:
	DeBase();
	
	void generateTrialSolution(const std::vector<Individual>& pop, size_t targetIndex, Individual& trial);
	bool selectSchema(const std::string& s);
	const std::string& getSchema() const;
	
protected:
	// base colution selection methods
	void rand(size_t& baseIndex) const;
	
	// difference selection choices 
	void one(const std::vector<Individual>& pop, size_t baseIndex, Individual& dest) const;
	void two(const std::vector<Individual>& pop, size_t baseIndex, Individual& dest) const;
	
	// recombination methods
	void exp(const Individual& parent, Individual& dest) const;
	void bin(const Individual& parent, Individual& dest) const;
};


template<class Individual, class Rand>
DeBase<Individual, Rand>::DeBase() {
	selectBase = &DeBase::rand;
	generateCandidate = &DeBase::one;
	recombineCandidate = &DeBase::bin;
}


template<class Individual, class Rand>
void DeBase<Individual, Rand>::generateTrialSolution(const std::vector<Individual>& pop, size_t targetIndex, Individual& trial) {
	size_t baseIndex;
	(this->*this->selectBase)(baseIndex);
	(this->*this->generateCandidate)(pop, baseIndex, trial);
	(this->*this->recombineCandidate)(pop[targetIndex], trial);
	static std::ofstream prva("prva.txt", std::ios::app);
	for (size_t i = 0; i < pop[targetIndex].chromosome.size(); ++i) {
		if (!std::isfinite(pop[targetIndex].chromosome[i])) {
			prva << ".bin " << targetIndex << ", " << pop.size() << " " << pop[targetIndex].chromosome << std::endl;
			break;
		}
	}
	trial.violation = -1;
}


template<class Individual, class Rand>
bool DeBase<Individual, Rand>::selectSchema(const std::string& s) {
	std::istringstream ss(s);
	
	static std::map<std::string, SelectBase> selectBaseMap;
	static std::map<std::string, GenerateCandidate> generateCandidateMap;
	static std::map<std::string, RecombineCandidate> recombineCandidateMap;
	
	if (selectBaseMap.empty()) {
		selectBaseMap["rand"] = &DeBase::rand;
		
		generateCandidateMap["1"] = &DeBase::one;
		generateCandidateMap["2"] = &DeBase::two;
		
		recombineCandidateMap["bin"] = &DeBase::bin;
		recombineCandidateMap["exp"] = &DeBase::exp;
	}
	
	std::string temp;
	bool noError = true;
	{
		{
			std::getline(ss, temp, '/');
			typename std::map<std::string, SelectBase>::const_iterator it = selectBaseMap.find(temp);
			if (it != selectBaseMap.end()) {
				schema = it->first;
				selectBase = it->second;
			} else {
				schema = selectBaseMap.begin()->first;
				selectBase = selectBaseMap.begin()->second;
				noError = false;
			}
		}
		{
			std::getline(ss, temp, '/');
			typename std::map<std::string, GenerateCandidate>::const_iterator it = generateCandidateMap.find(temp);
			if (it != generateCandidateMap.end()) {
				schema += "/" + it->first;
				generateCandidate = it->second;
			} else {
				schema += "/" + generateCandidateMap.begin()->first;
				generateCandidate = generateCandidateMap.begin()->second;
				noError = false;
			}
		}
		{
			std::getline(ss, temp);
			typename std::map<std::string, RecombineCandidate>::const_iterator it = recombineCandidateMap.find(temp);
			if (it != recombineCandidateMap.end()) {
				schema += "/" + it->first;
				recombineCandidate = it->second;
			} else {
				schema += "/" + recombineCandidateMap.begin()->first;
				recombineCandidate = recombineCandidateMap.begin()->second;
				noError = false;
			}
		}
	}
	
	return noError;
}


template<class Individual, class Rand>
const std::string& DeBase<Individual, Rand>::getSchema() const {
	return schema;
}
	

template<class Individual, class Rand>
void DeBase<Individual, Rand>::rand(size_t& baseIndex) const {
	baseIndex = Rand().exclusiveInterval(indices.size());
}


template<class Individual, class Rand>
void DeBase<Individual, Rand>::one(const std::vector<Individual>& pop, size_t baseIndex, Individual& dest) const {
	assert(scalingFactor.size() >= 1);
	
	// random parent and vector difference selection (result is 3 mutually different indices)
	size_t p1 = baseIndex;
	size_t p2 = p1;
	while (p2 == p1) 
		p2 = Rand().exclusiveInterval(pop.size());
	size_t p3 = p1;
	while ((p2 == p3) || (p1 == p3))
		p3 = Rand().exclusiveInterval(pop.size());
	
	dest.chromosome = pop[p1].chromosome + (pop[p2].chromosome - pop[p3].chromosome) * scalingFactor[0];
}


template<class Individual, class Rand>
void DeBase<Individual, Rand>::two(const std::vector<Individual>& pop, size_t baseIndex, Individual& dest) const {
	assert(pop.size() > 5);
	assert(scalingFactor.size() >= 2);
	
	// random parent and vector difference selection (result is 5 mutually different indices)
	size_t p1 = baseIndex;
	size_t p2 = p1;
	while (p2 == p1) 
		p2 = Rand().exclusiveInterval(pop.size());
	size_t p3 = p1;
	while ((p2 == p3) || (p1 == p3))
		p3 = Rand().exclusiveInterval(pop.size());
	size_t p4 = p1;
	while ((p4 == p3) || (p4 == p2) || (p4 == p1)) 
		p4 = Rand().exclusiveInterval(pop.size());
	size_t p5 = p1;
	while ((p5 == p4) || (p5 == p3) || (p5 == p2) || (p5 == p1)) 
		p5 = Rand().exclusiveInterval(pop.size());
	
	dest.chromosome = pop[p1].chromosome + (pop[p2].chromosome - pop[p3].chromosome) * scalingFactor[0] + 
		(pop[p4].chromosome - pop[p5].chromosome) * scalingFactor[1];
}


template<class Individual, class Rand>
void DeBase<Individual, Rand>::exp(const Individual& parent, Individual& dest) const {
	// trial genes are set as follows:
	// random gene is taken from mutant (which is provided as an ininital value to 'dest')
	// going from that randomly selected gene in circular path around all genes,
	//		genes are taken from mutant while rand(0,1) < Cr, after first random fails
	//		the test, all remaining genes are taken from parent
	size_t start = Rand().exclusiveInterval(dest.chromosome.size());
	bool cross = true;
	for (size_t i = start + 1; (i != start) && (i != start + dest.chromosome.size()); ++i) {
		if (i >= dest.chromosome.size())
			i = 0;
		
		if (!cross || (Rand().interval(1.0) > crossoverProbability)) {
			cross = false;
			dest.chromosome[i] = parent.chromosome[i];
		}
	}
}


template<class Individual, class Rand>
void DeBase<Individual, Rand>::bin(const Individual& parent, Individual& dest) const {
	size_t start = Rand().exclusiveInterval(dest.chromosome.size());
	for (size_t i = 0; i < dest.chromosome.size(); ++i) {
		if ((i!=start) && (Rand().interval(1.0) > crossoverProbability)) 
			dest.chromosome[i] = parent.chromosome[i];
	}
}


#endif // DEBASE_H_INCLUDED
