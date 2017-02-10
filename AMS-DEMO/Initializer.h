#ifndef INITIALIZER_H_INCLUDED
#define INITIALIZER_H_INCLUDED


#include <Random.h>
#include <cassert>
#include "TypeWrapper.h"


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// template class DynamicRandomInitializer
/// used as a parameter for NumericOptimizer<...>::init(initializer);
///
/// uses Random::CRand to generate random solution vectors
/// CRand seed should be set before first use
template<class Float, class Solution>
class DynamicRandomInitializer {
	std::vector<Float> min;
	std::vector<Float> max;
	size_t chromosomeSize;
	size_t numOfCriteria;
	size_t numOfProperties;
	mutable size_t preLoadedUsed;
	
public:
	std::vector<Solution> preLoaded;
	
public:
	DynamicRandomInitializer(size_t chSize, size_t critSize, size_t numProp = 0);
	
	void generateSolution(Solution& sol) const;
	void defineGeneValueRange(size_t geneIndex, double vMin, double vMax);
};


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// template class DynamicRandomInitializer
template<class Float, class Solution>
DynamicRandomInitializer<Float, Solution>::DynamicRandomInitializer(size_t chSize, size_t critSize, size_t numProp) {
	chromosomeSize = chSize;
	numOfCriteria = critSize;
	numOfProperties = numProp;
	min.resize(chromosomeSize, 0);
	max.resize(chromosomeSize, 0);
	preLoadedUsed = 0;
}


template<class Float, class Solution>
void DynamicRandomInitializer<Float, Solution>::generateSolution(Solution& sol) const {
	if (TypeWrapper::size(sol.chromosome) != chromosomeSize) {
		TypeWrapper::resize(sol.chromosome, chromosomeSize);
	}
	if (TypeWrapper::size(sol.criteria) != numOfCriteria) {
		TypeWrapper::resize(sol.criteria, numOfCriteria);
	}
	if (sol.hasProperties && (TypeWrapper::size(sol.getProperties()) != numOfProperties)) {
		TypeWrapper::resize(sol.criteria, numOfProperties);
	}
	
	if ((preLoadedUsed < preLoaded.size()) && TypeWrapper::size(preLoaded[preLoadedUsed].chromosome) == chromosomeSize) {
		sol = preLoaded[preLoadedUsed];
		++preLoadedUsed;
	} else {	
		for (size_t i = 0; i < sol.chromosome.size(); ++i) {
			assert(min[i] <= max[i]);
			sol.chromosome[i] = Random::CRand().interval(min[i], max[i]);
			sol.violation = -1;
		}
	}
}


template<class Float, class Solution>
void DynamicRandomInitializer<Float, Solution>::defineGeneValueRange(size_t geneIndex, double vMin, double vMax) {
	assert(geneIndex < chromosomeSize);
	
	min[geneIndex] = vMin;
	max[geneIndex] = vMax;
}


#endif // INITIALIZER_H_INCLUDED
