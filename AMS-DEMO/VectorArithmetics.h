#ifndef VECTORARITHMETICS_H_INCLUDED
#define VECTORARITHMETICS_H_INCLUDED


#include <vector>
#include <cassert>
#include <iostream>


template<class C, class T>
std::basic_ostream<C>& operator<< (std::basic_ostream<C>& out, const std::vector<T>& op) {
	out << "<" ;
	if (op.size() > 0)
		out << op[0];
	for (size_t i = 1; i < op.size(); ++i)
		out << "," << op[i];
	out << ">";
	return out;
}


template<class T>
std::vector<T> operator- (const std::vector<T>& a, const std::vector<T>& b) {
	assert(a.size() == b.size());
	std::vector<T> c(a.size());
	for (size_t i = 0; i < c.size(); ++i) {
		c[i] = a[i] - b[i];
	}
	return c;
}


template<class T>
std::vector<T> operator+ (const std::vector<T>& a, const std::vector<T>& b) {
	assert(a.size() == b.size());
	std::vector<T> c(a.size());
	for (size_t i = 0; i < c.size(); ++i) {
		c[i] = a[i] + b[i];
	}
	return c;
}


template<class T, class S>
std::vector<T> operator* (const std::vector<T>& a, S scalar) {
	std::vector<T> c(a.size());
	for (size_t i = 0; i < c.size(); ++i) {
		c[i] = a[i] * scalar;
	}
	return c;
}


#endif // VECTORARITHMETICS_H_INCLUDED
