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
