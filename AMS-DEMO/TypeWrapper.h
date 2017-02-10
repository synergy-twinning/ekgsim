#ifndef TYPEWRAPPER_H_INCLUDED
#define TYPEWRAPPER_H_INCLUDED


#include <vector>
#include <stdexcept>
#include "Array.h"


struct TypeWrapper {
	// char
	static size_t size(const char*) {return 0;}
	static void resize(const char*) {throw std::runtime_error("Error: resize called on character");}

	// vector
	template<class T>
	static size_t size(const std::vector<T>& a) {return a.size();}

	template<class T>
	static void resize(std::vector<T>& a, size_t n) {a.resize(n);}
	
	// array
	template<class T, size_t N>
	static size_t size(const Array<T, N>& a) {return a.size();}

	template<class T, size_t N>
	static void resize(const Array<T, N>& a, size_t n) {
		if (n != N)
			throw std::runtime_error("error: trying to resize Array<T, N>");
	}
};


#endif // TYPEWRAPPER_H_INCLUDED
