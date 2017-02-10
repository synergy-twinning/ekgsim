#ifndef DEMO_H_INCLUDED
#define DEMO_H_INCLUDED


#include "De.h"


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// template class Demo
/// used as an Algorithm policy for NumericOptimizer<...>;
///
template<class Chromosome, class Evaluation>
class Demo : 
	// inherit base DE definitions
	public DeBase<IndividualStruc<Chromosome, typename Evaluation::Value> >, 
	// make De conform the form of NumericOptimizer
	public NumericOptimizer<Demo<Chromosome, Evaluation>, Chromosome, Evaluation> 
{
protected:
	void (Demo::*truncatePopulation)(std::vector<size_t>& indices);
	friend class NumericOptimizer<Demo<Chromosome, Evaluation>, Chromosome, Evaluation>;
	
public:
	// algorithm parameters
	bool alwaysAdd; // no selection is performed until the next generation
	
public:
	Demo();
	
	// returnes indices of individuals that form pareto front
	void getFront(std::vector<size_t>& front);
	
protected:
	void nextGeneration();
	void generationalSelection();
	
	void performSelection(size_t targetIndex, size_t trialIndex);
	void truncatePopulationSpea2(std::vector<size_t>& indices);
};


template<class Chromosome, class Evaluation>
Demo<Chromosome, Evaluation>::Demo() {
	alwaysAdd = false;
	truncatePopulation = &Demo::truncatePopulationSpea2;
}


template<class Chromosome, class Evaluation>
void Demo<Chromosome, Evaluation>::getFront(std::vector<size_t>& front) {
	size_t numInd = this->pop.size();
	std::vector<size_t> dominatedCount(numInd, 0);
	
	for (size_t i = 0; i < numInd; ++i) {
		for (size_t j = 0; j < numInd; ++j) {
			if (this->pop[i] < this->pop[j])
				++dominatedCount[j];
		}
	}
	front.reserve(numInd);
	front.clear();
	for (size_t i = 0; i < numInd; ++i) {
		if (dominatedCount[i] == 0) {
			front.push_back(i);
		}
	}
}


template<class Chromosome, class Evaluation>
void Demo<Chromosome, Evaluation>::nextGeneration() {
	randomizeIndices(this->indices, this->pop.size());
	
	for (size_t i = 0; i < this->indices.size(); ++i) {
		this->generateTrialSolution(this->pop, this->indices[i], this->children[i]);
		this->eval.normalize(this->children[i].chromosome);
	}
}


template<class Chromosome, class Evaluation>
void Demo<Chromosome, Evaluation>::generationalSelection() {
	// create joint population "pop" from "pop" and "children"
	for (size_t i = 0; i < this->indices.size(); ++i) {
		performSelection(this->indices[i], i);
	}
	
	// truncate "pop" -> produce indices of surviving individuals of population "pop"
	(this->*truncatePopulation)(this->indices);
	
	// use population "children" for temporary storage
	for (size_t i = 0; i < this->indices.size(); ++i) {
		this->children[i] = this->pop[this->indices[i]];
	}
	
	// put resulting population and its objective values back to "pop" and "popVal"
	this->pop.swap(this->children);
	this->children.resize(this->indices.size());
}


template<class Chromosome, class Evaluation>
void Demo<Chromosome, Evaluation>::performSelection(size_t targetIndex, size_t trialIndex) {
	if (!alwaysAdd && (this->children[trialIndex] <= this->pop[targetIndex])) {
		// child dominates parent
		this->pop[targetIndex] = this->children[trialIndex];
	} else if (!alwaysAdd && (this->pop[targetIndex] <= this->children[trialIndex])) {
		// parent dominates child
	} else {
		// neither parent or child dominates the other
		this->pop.push_back(this->children[trialIndex]);
	}
}


template<class Chromosome, class Evaluation>
void Demo<Chromosome, Evaluation>::truncatePopulationSpea2(std::vector<size_t>& indices) {
	size_t numberOfIndividuals = this->pop.size();
	size_t popSize = indices.size();
	int k = (int)sqrt(numberOfIndividuals);
	//std::cout << "pop size = " << numberOfIndividuals << " / " << popSize << " \n";
	
	if (popSize == numberOfIndividuals) {
		// short circuit for border case scenario of no extra individuals
		for (size_t i = 0; i < indices.size(); ++i)
			indices[i] = i;
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
			indices.resize(last + 1);
			for (size_t i = 0; i <= last; ++i)
				indices[i] = sortedFitnessI[i];
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
			if (numRedundant > 0) std::cerr << "truncating\n";
			std::vector<size_t> marked;
			for (size_t j = 0; j < numRedundant; ++j) {
				// put candidates for deletition into vector marked
				marked.clear();
				marked.resize(indices.size(), 0);
				
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
				std::cerr << "\n";
				for (size_t c = 0; c < count; ++c)
					std::cerr << marked[c] << " "; 
				std::cerr << "\n";
				
				size_t deleteIndex = marked[Random::CRand().exclusiveInterval(count)];
				// remove from record
				for (size_t jj = 0; jj < numberOfIndividuals; ++jj) {
					if ((distances[deleteIndex][jj] == 0) && (copies[jj] > 0))
						--copies[jj];
				}
				copies[deleteIndex] = 0;
				
				assert(std::find(indices.begin(), indices.end(), deleteIndex) != indices.end());
				indices.erase(std::find(indices.begin(), indices.end(), deleteIndex));
				
				std::cerr << " " << this->pop[deleteIndex].criteria;
				double distAvg = 0.0, distAvgD = 0.0;
				size_t cnt = 0, cntD = 0;
				for (size_t i = 0; i < numberOfIndividuals; ++i) {
					distAvg += distances[i][k];
					cnt++;
				}
			}
				
		} else {
			indices.swap(sortedFitnessI);
			indices.resize(popSize);
		}
	}
}


#endif // DEMO_H_INCLUDED
