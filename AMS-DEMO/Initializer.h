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
