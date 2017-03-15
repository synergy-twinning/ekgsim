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

#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED


#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <Random.h>
#include "timer.h"


void randomizeIndices(std::vector<size_t>& indices, size_t size = 0);


template<class T>
struct SortItem {
	size_t i;
	T val;
	
	inline operator size_t() const {
		return i;
	}
	
	friend inline bool operator< (const SortItem& a, const SortItem& b) {
		return a.val < b.val;
	}
};


template<class T>
void getSortedIndices(const std::vector<T>& src, std::vector<size_t>& indices) {
	std::vector<SortItem<T> > items(src.size());
	
	for (size_t i = 0; i < src.size(); ++i) {
		items[i].val = src[i];
		items[i].i = i;
	}
	
	std::sort(items.begin(), items.end());
	
	indices.resize(items.size());
	std::copy(items.begin(), items.end(), indices.begin());
}


template<class T>
bool areIndicesOk(const std::vector<T>& indices, size_t maxIndex) {
	std::vector<T> copy = indices;
	std::sort(copy.begin(), copy.end());
	
	bool ok = (copy[0] >= 0);
	if (!ok) std::cout << "sorted_indices[0]=" << copy[0] << " ";
	bool ok1 = (copy[copy.size()-1] < maxIndex);
	if (!ok1) std::cout << "sorted_indices[last]=" << copy[copy.size()-1] << " ";
	ok &= ok1;
	for (size_t i = 1; i < copy.size(); ++i) {
		bool ok2 = (copy[i] != copy[i-1]);
		if (!ok2) std::cout << "sorted_indices " << (i-1) << " & " << i << " = " << copy[i] << " ";
		ok &= ok2;
	}
		
	if (!ok) std::cout << "\n";
	return ok;
}


#endif // UTILITIES_H_INCLUDED
