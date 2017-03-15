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

#ifndef PDEMO_H_INCLUDED
#define PDEMO_H_INCLUDED


#include "ParallelNumericOptimizer.h"
#include "DeBase.h"


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// template class PDemo
/// used as an Algorithm policy for NumericOptimizer<...>;
///
template<class Chromosome, class Evaluation>
class PDemo : 
	// inherit base DE definitions
	public DeBase<IndividualStruc<Chromosome, typename Evaluation::Value, typename Evaluation::Properties> >, 
	// make PDemo conform the form of NumericOptimizer
	public ParallelNumericOptimizer<PDemo<Chromosome, Evaluation>, Chromosome, Evaluation> 
{
protected:
	void (PDemo::*truncatePopulation)(std::vector<size_t>& indices);
	friend class ParallelNumericOptimizer<PDemo<Chromosome, Evaluation>, Chromosome, Evaluation>;
	size_t parentIndex; // used in generateCandidate function
	
public:
	typedef IndividualStruc<Chromosome, typename Evaluation::Value, typename Evaluation::Properties> Individual;
	bool alwaysAdd;
	
public:
	PDemo();
	
	// returnes indices of individuals that form pareto front
	void getFront(std::vector<size_t>& front) const;
	
protected:
	void nextGeneration();
	void performTruncate();
	
	void performSelection(size_t targetIndex, Individual trial);
	void truncatePopulationSpea2(std::vector<size_t>& indices);
	
	Chromosome generateCandidate();
};


template<class Chromosome, class Evaluation>
PDemo<Chromosome, Evaluation>::PDemo() : parentIndex(0) {
	alwaysAdd = false;
	truncatePopulation = &PDemo::truncatePopulationSpea2;
}


template<class Chromosome, class Evaluation>
void PDemo<Chromosome, Evaluation>::getFront(std::vector<size_t>& front) const {
	size_t numInd = this->pop.size();
	std::vector<size_t> dominatedCount(numInd, 0);
	
	for (size_t i = 0; i < numInd; ++i) {
		for (size_t j = 0; j < numInd; ++j) {
			if (this->pop[i] < this->pop[j])
				++dominatedCount[j];
		}
	}
	front.clear();
	front.reserve(numInd);
	for (size_t i = 0; i < numInd; ++i) {
		if (dominatedCount[i] == 0) {
			front.push_back(i);
		}
	}
}


template<class Chromosome, class Evaluation>
void PDemo<Chromosome, Evaluation>::performTruncate() {
	// truncate "pop" -> produce indices of surviving individuals of population "pop"
	randomizeIndices(this->indices, this->populationSize);
	(this->*truncatePopulation)(this->indices);
	
	// use temporary storage population with organised individuals
	std::vector<Individual> newPop(this->indices.size());
	for (size_t i = 0; i < this->indices.size(); ++i) {
		newPop[i] = this->pop[this->indices[i]];
	}
	
	// put resulting population and its objective values back to "pop" 
	this->pop.swap(newPop);
}


template<class Chromosome, class Evaluation>
void PDemo<Chromosome, Evaluation>::performSelection(size_t targetIndex, Individual trial) {
	if (!alwaysAdd && (trial <= this->pop[targetIndex])) {
		// child dominates parent
		this->pop[targetIndex] = trial;
	} else if (!alwaysAdd && (this->pop[targetIndex] <= trial)) {
		// parent dominates child
		// do nothing
	} else {
		// neither parent or child dominates the other
		this->pop.push_back(trial);
	}
}


template<class Chromosome, class Evaluation>
void PDemo<Chromosome, Evaluation>::truncatePopulationSpea2(std::vector<size_t>& truncIndices) {
	size_t numberOfIndividuals = this->pop.size();
	size_t popSize = truncIndices.size();
	int k = (int)sqrt((double)numberOfIndividuals);
	//std::cout << "pop size = " << numberOfIndividuals << " / " << popSize << " \n";
	
	if (popSize == numberOfIndividuals) {
		// short circuit for border case scenario of no extra individuals
		for (size_t i = 0; i < truncIndices.size(); ++i)
			truncIndices[i] = i;
	} else {
		// calculate raw fitness
		std::vector<double> rawFitness(numberOfIndividuals, 0);
		{
			std::vector<size_t> strength(numberOfIndividuals, 0);
			for (size_t i = 0; i < numberOfIndividuals; ++i) {
				for (size_t j = 0; j < numberOfIndividuals; ++j) {
					if (this->pop[i] < this->pop[j])
						++strength[i];
				}
			}
			
			for (size_t i = 0; i < numberOfIndividuals; ++i) {
				for (size_t j = 0; j < numberOfIndividuals; ++j) {
					if (this->pop[j] < this->pop[i])
						rawFitness[i] += strength[j];
				}
			}
		}
		
		// calculate distance
		std::vector<std::vector<double> > distances;
		std::vector<size_t> copies(numberOfIndividuals, 1);
		std::vector<std::vector<int> > nn(numberOfIndividuals);
		
		for (size_t i = 0; i < numberOfIndividuals; ++i)
			nn[i].resize(numberOfIndividuals, -1);
		
		std::vector<double> tmp(numberOfIndividuals, 0);
		for (size_t i = 0; i < numberOfIndividuals; ++i) {
			nn[i][0] = i;
			for (size_t j = 0; j < i; ++j) {
				assert(distances.size() > j);
				assert(distances[j].size() > i);
				tmp[j] = distances[j][i];
			}
			tmp[i] = 0;
			
			for (size_t j = i + 1; j < numberOfIndividuals; ++j) {
				tmp[j] = vectorDistance(this->pop[i].criteria, this->pop[j].criteria);
				
				if (tmp[j] == 0) {
					nn[i][copies[i]] = j;
					nn[j][copies[j]] = i;
					++copies[i];
					++copies[j];
				}
			}
			distances.push_back(tmp);
		}
		
		// calculate fitness
		std::vector<double> fitness(numberOfIndividuals);
		for (size_t i = 0; i < numberOfIndividuals; ++i) {
			fitness[i] = rawFitness[i] + 1.0 / (2.0 + distances[i][k]);
		}
		
		// get sorted indices of fitness vector
		std::vector<size_t> sortedFitnessI;
		getSortedIndices(fitness, sortedFitnessI);
		
		if (fitness[sortedFitnessI[popSize]] < 1.0) {
			
			size_t last = popSize;
			for (; (last < numberOfIndividuals-1) && (fitness[sortedFitnessI[last + 1]] < 1.0); ++last);
			truncIndices.resize(last + 1);
			for (size_t i = 0; i <= last; ++i)
				truncIndices[i] = sortedFitnessI[i];
			for (size_t i = last + 1; i < numberOfIndividuals; ++i) {
				// remove *from record*
				for (size_t j = 0; j < numberOfIndividuals; ++j) {
					if ((distances[sortedFitnessI[i]][j] == 0) && (copies[j] > 0))
						--copies[j];
				}
				copies[sortedFitnessI[i]] = 0;
			}
			
			// truncate population
			size_t numRedundant = last - popSize + 1;
//			if (numRedundant > 0) std::cerr << "truncating\n";
			std::vector<size_t> marked;
			for (size_t j = 0; j < numRedundant; ++j) {
				// put candidates for deletition into vector marked
				marked.clear();
				marked.resize(truncIndices.size(), 0);
				
				// ? first, mark individuals that have copies (other individuals with 
				// the same score in criteria); put down the ones with most copies first
				//   count contains ~ the number of individuals with copies
				//   maxCopies contains the max number of identical individuals
				size_t maxCopies = 0;
				size_t count = 0;
				for (size_t i = 0; i < numberOfIndividuals; ++i) {
					if (copies[i] > maxCopies) {
						count = 0;
						maxCopies = copies[i];
					}
					if (copies[i] == maxCopies) {
						marked[count] = i;
						++count;
					}
				}
				
				// ? there is more individuals with the same number of copies than there is 
				// copies of single individual
				// (usually this happens when there is no copies at all; count then holds 
				// total number of individuals nad maxCopies = 1)
				if (count > maxCopies) {
					std::vector<size_t> neighbour(count, 1);
					for (; count > maxCopies;) {
						double minDistance = 1e32;
						size_t count2 = 0;
						
						for (size_t i = 0; i < count; ++i) {
							double distance = -1.0;
							
							while ((distance < 0.0) && (neighbour[i] < numberOfIndividuals)) {
								// get distance to NN (between marked[i] and neighbour[i])
								// int neighborIndex = GetNN(marked[i], neighbour[i]);
								
								if (nn[marked[i]][neighbour[i]] == -1) {
									double minDistance1 = 2e32;
									int minIndex = 0;
									size_t prevMinIndex = nn[marked[i]][neighbour[i] - 1];
									double prevMinDistance = distances[marked[i]][prevMinIndex];
									for (size_t ii = 0; ii < numberOfIndividuals; ii++) {
										double distance1 = distances[marked[i]][ii];
										
										if ((distance1 < minDistance1) && (marked[i] != ii)) {
											if ((distance1 > prevMinDistance) || 
												((distance1 == prevMinDistance) && (ii > prevMinIndex))) 
											{
												minDistance1 = distance1;
												minIndex = ii;
											}
										}
									}
									
									nn[marked[i]][neighbour[i]] = minIndex;
								}
								int neighbourIndex = nn[marked[i]][neighbour[i]];
								if (copies[neighbourIndex] != 0)
									distance = distances[marked[i]][neighbourIndex];
								else
									distance = -1;
								
								++neighbour[i];
							}
							
							if (distance < minDistance) {
								count2 = 0;
								minDistance = distance;
							}
							if (distance == minDistance) {
								marked[count2] = marked[i];
								neighbour[count2] = neighbour[i];
								++count2;
							}
						}
						
						count = count2;
						if (minDistance == -1)
							break;
					}
				}
				
				size_t deleteIndex = marked[Random::CRand().exclusiveInterval(count)];
				// remove from record
				for (size_t jj = 0; jj < numberOfIndividuals; ++jj) {
					if ((distances[deleteIndex][jj] == 0) && (copies[jj] > 0))
						--copies[jj];
				}
				copies[deleteIndex] = 0;
				
				assert(std::find(truncIndices.begin(), truncIndices.end(), deleteIndex) != truncIndices.end());
				truncIndices.erase(std::find(truncIndices.begin(), truncIndices.end(), deleteIndex));
				
				double distAvg = 0.0;//, distAvgD = 0.0;
				size_t cnt = 0;//, cntD = 0;
				for (size_t i = 0; i < numberOfIndividuals; ++i) {
					distAvg += distances[i][k];
					cnt++;
				}
			}
				
		} else {
			truncIndices.swap(sortedFitnessI);
			truncIndices.resize(popSize);
		}
	}
}


template<class Chromosome, class Evaluation>
Chromosome PDemo<Chromosome, Evaluation>::generateCandidate() {
	Individual ret;
	if (parentIndex >= this->pop.size())
		parentIndex = 0;
	this->generateTrialSolution(this->pop, parentIndex, ret);
	++parentIndex;
	return ret.chromosome;
};

#endif // PDEMO_H_INCLUDED
