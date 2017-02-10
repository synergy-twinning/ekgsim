#include "utilities.h"
#include <iostream>


struct SortElement {
	size_t index;
	size_t rand;
	
	friend bool operator< (const SortElement& a, const SortElement& b) {return a.rand < b.rand;}
	operator size_t() const {return index;}
};


void randomizeIndices(std::vector<size_t>& indices, size_t size) {
	// indices might not be defined yet
	if (size != 0) 
		indices.resize(size);
	
	std::vector<SortElement> randoms(indices.size());
	for (size_t i = 0; i < randoms.size(); ++i) {
		randoms[i].index = i;
		randoms[i].rand = Random::CRand();
	}
	
	std::sort(randoms.begin(), randoms.end());
	std::copy(randoms.begin(), randoms.end(), indices.begin());
}
