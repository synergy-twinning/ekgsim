#ifndef HYPERVOLUME_H_INCLUDED
#define HYPERVOLUME_H_INCLUDED


#include <vector>
#include <algorithm>
#include "utilities.h"
#include <iostream>
#include <stdexcept>


template<class Individual>
struct HsoSlice {
	// 1-dimensional volume of the current slice. 
	double depth;
	std::vector<Individual> individuals;

	HsoSlice(double d, const std::vector<Individual>& ind) : depth(d) {
		individuals = ind;
	}
	
	const HsoSlice& operator= (const HsoSlice& slice) {
		depth = slice.depth;
		individuals = slice.individuals;
		return *this;
	}
};


template<class Individual>
class HsoPopulation {
protected:
	// non-dominated individuals
	std::vector<Individual> individuals;
	std::vector<double> referencePoint;
	
public:
	const std::vector<double>& criteriaMin;
	const std::vector<double>& criteriaMax;
	
	HsoPopulation(const std::vector<Individual>& pop, const std::vector<double>& cmin, const std::vector<double>& cmax) :
		criteriaMin(cmin), criteriaMax(cmax)
	{
		std::vector<size_t> n(pop.size(), 0);
		std::vector<size_t> eliminate;
		
		for (size_t i = 0; i < pop.size(); i++) {
			for (size_t j = 0; j < pop.size(); j++) {
				if ((i != j) && (pop[j].violation == 0)) {
					if (pop[j].criteria < pop[i].criteria)
						n[i]++;
				}
			}
			if ((n[i] == 0) && (pop[i].violation == 0)) {
				individuals.push_back(pop[i]);
			}
		}
		
		referencePoint.resize(individuals[0].criteria.size(), 2);
	}
	
	double getHypMeasure() {
		if (individuals.size() == 0) {
			return 0;
		}
		
		scaleCriteria();
		sortObjective(0);
		
		size_t n = individuals[0].criteria.size();
		std::vector<HsoSlice<Individual> > slices;
		std::vector<HsoSlice<Individual> > newSlices;
		std::vector<HsoSlice<Individual> > nextSlices;
		slices.push_back(HsoSlice<Individual>(1, individuals));
		
		for (size_t j = 0; j < n - 1; j++) {
			newSlices.clear();
			for (size_t i = 0; i < slices.size(); i++) {
				getSlices(slices[i].individuals, j, nextSlices);
				for (size_t k = 0; k < nextSlices.size(); k++) {
					newSlices.push_back(HsoSlice<Individual>(slices[i].depth * nextSlices[k].depth,
						nextSlices[k].individuals));
				}
			}
			slices.swap(newSlices);
		}
		
		double volume = 0;
		for (size_t i = 0; i < slices.size(); i++) {
            /*
			volume += slices[i].depth * fabs(std::min(0.0, slices[i].individuals[0].criteria[n - 1] - 
				referencePoint[n - 1]));
                * */
            volume += slices[i].depth * 
                std::max(0.0, referencePoint[n - 1] - slices[i].individuals[0].criteria[n - 1]);
		}
		
		return volume;
	}

	void scaleCriteria() {
		if ((individuals[0].criteria.size() > criteriaMin.size()) || 
			(individuals[0].criteria.size() > criteriaMax.size())) 
		{
			throw std::runtime_error("The number of criteria of individual does not match the number of predefined\
				criteria minimas or maximas for hypervolume");
		}
		
		for (size_t i = 0; i < individuals.size(); i++) {
			for (size_t j = 0; j < individuals[0].criteria.size(); j++) {
                if (individuals[i].criteria[j] > criteriaMax[j])
                    individuals[i].criteria[j] = criteriaMax[j];
                individuals[i].criteria[j] = (individuals[i].criteria[j] - criteriaMin[j]) /
                        (criteriaMax[j] - criteriaMin[j]);
			}
		}
	}

	void sortObjective(size_t index) {
		std::vector<size_t> sortIndex(individuals.size());
		std::vector<double> sortValue(individuals.size());
		std::vector<Individual> sortedInd(individuals.size());
		
		for (size_t i = 0; i < individuals.size(); i++)
			sortValue[i] = individuals[i].criteria[index];	
		getSortedIndices(sortValue, sortIndex);
		
		for (size_t i = 0; i < individuals.size(); i++)
			sortedInd[i] = individuals[sortIndex[i]];
		individuals.swap(sortedInd);
	}

	void getSlices(std::vector<Individual>& inds, size_t index, std::vector<HsoSlice<Individual> >& slice) {
		if (inds.size() < 1)
			throw std::runtime_error("error in HsoPop.getSlices: empty vector inds");
		
		std::vector<Individual> newInds;
		typename std::vector<Individual>::iterator point = inds.begin();
		slice.clear();
		
		typename std::vector<Individual>::iterator newPoint = point + 1;
		for (; newPoint != inds.end(); newPoint = point + 1) {
			newInds = insert(*point, index + 1, newInds);
			slice.push_back(HsoSlice<Individual>(fabs(point->criteria[index] - 
				newPoint->criteria[index]), newInds));
			point = newPoint;
		};
		
		newInds = insert(*point, index + 1, newInds);
		slice.push_back(HsoSlice<Individual>(fabs(point->criteria[index] - 
			referencePoint[index]), newInds));
	}

	std::vector<Individual> insert(const Individual& point, size_t depth, const std::vector<Individual>& inds) {
		std::vector<Individual> result;
		
		size_t i = 0;
		for (; (i < inds.size()) && (inds[i].criteria[depth] < point.criteria[depth]); ++i) {
			result.push_back(inds[i]);
		}
		result.push_back(point);
		
		for (; i < inds.size(); ++i) {
			if (!dominates(point, inds[i], depth)) {
				result.push_back(inds[i]);
			}
		}
		
		return result;
	}
	
	bool dominates(const Individual& a, const Individual& b, size_t k) {
		bool dom 	= (a.criteria[0] < b.criteria[0]);
		bool cover 	= (a.criteria[0] <= b.criteria[0]);
		for (; k < a.criteria.size(); ++k) {
			dom 	|= (a.criteria[k] < b.criteria[k]);
			cover 	&= (a.criteria[k] <= b.criteria[k]);
		}
		return dom & cover;
	}
};

#endif // HYPERVOLUME_H_INCLUDED
