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

#ifndef VECTORMATH_H_INCLUDED
#define VECTORMATH_H_INCLUDED

/** @file vectorMath.h
    @brief routines for vector arithmetics
**/

#include <vector>
#include <cmath>
#include <algorithm>

#include <iostream>
#define DBVAL(A) std::cout << #A << " = " << (A) << " ";
#define DBVALn(A) DBVAL(A); std::cout << "\n";


typedef std::vector<double> TFunction;


/// returns sign of the scalar number (either 1 or -1)
template<class T>
T sgn(T num) {
	return (num >= 0 ? T(1) : T(-1));
}


/// returns the vector containing squares of the elements (in a range from start to end, end not included) of the original vector (matlab: src .^ 2)
template<class T>
std::vector<T> sqr(const std::vector<T>& src, size_t start, size_t end) {
	std::vector<T> s(src.begin() + start, src.begin() + end);
	for (size_t i = 0; i < end - start; ++i) {
		s[i] *= s[i];
	}
	return s;
}


/// returns the vector containing squares of the elements of the original vector (matlab: src .^ 2)
template<class T>
std::vector<T> sqr(const std::vector<T>& src) {
	return sqr(src, 0, src.size());
}


/// inline multiplication of vector src with scalar m
template<class T, class S>
void mult(std::vector<T>& src, S m) {
	for (typename std::vector<T>::iterator it = src.begin(); it != src.end(); ++it)
		*it *= m;
}


/// inline addition of vector src and scalar a
template<class T, class S>
void add(std::vector<T>& src, S a) {
	for (typename std::vector<T>::iterator it = src.begin(); it != src.end(); ++it)
		*it += a;
}


/// derives vector src, saves result in vector dest
template<class T>
void derive(const std::vector<T>& src, std::vector<T>& dest) {
	dest.resize(src.size()-1);
	for (size_t i = 0; i < dest.size(); ++i) {
		dest[i] = src[i+1] - src[i];
	}
}


/// returns indices to elements preceding the sign change
template<class T, class Int>
void findZeros(const std::vector<T>& src, std::vector<Int>& zeros) {
	T prevSgn = sgn(src[0]);
	for (size_t i = 1; i < src.size(); ++i) {
		T newSgn = sgn(src[i]);
		if (prevSgn != newSgn) {
			prevSgn = newSgn;
			zeros.push_back(i);
		}
	}
}


/// dot product of two vectors
template<class T>
T dot(const std::vector<T>& a, size_t startA, const std::vector<T>& b, size_t startB, size_t len) {
	T s = 0;
	for (size_t i = startA, j = startB; i < (startA + len); ++i, ++j) {
		s += a[i] * b[j];
	}
	return s;
}


/// summation of the vector elements (in the range from start to exclusive end)
template<class T>
T sum(const std::vector<T>& src, size_t start, size_t end) {
	T s = 0;
	for (size_t i = start; i < end; ++i) {
		s += src[i];
	}
	return s;
}


/// summation of the vector elements
template<class T>
T sum(const std::vector<T>& src) {
	return sum(src, 0, src.size());
}


/// calculats mean and variation of vector
template<class T>
void meanAndVar(T& mean, T& var, const std::vector<T>& src, size_t start, size_t end) {
	T inverseSize = T(1) / T(end - start);
	mean = sum(src, start, end) * inverseSize;
	var = 0;
	for (size_t i = start; i < end; ++i) {
		var += sqr(src[i] - mean);
	}
	var *= inverseSize;
}


/// calculates minimum and maximum of a vector
template<class T>
void minAndMax(T& minOut, T& maxOut, const std::vector<T>& src) {
	if (src.size() == 0)
		return;
	else
		minOut = maxOut = src[0];
		
	for (size_t i = 1; i < src.size(); ++i) {
		if (src[i] < minOut)
			minOut = src[i];
		else if (src[i] > maxOut)
			maxOut = src[i];
	}
}


/// linear resample (sampling frequency is multiplied by factor)
template<class T>
void resample(const std::vector<T>& src, std::vector<T>& dest, double factor) {
	if (factor == 1) {
		dest = src;
	} else {
		dest.resize(1u + (size_t)ceil((src.size() - 1) * factor));
		
		for (size_t i = 0; i < dest.size(); ++i) {
			double oldPos = i/factor;
			double fPos = floor(oldPos);
			double cPos = 1.0 + fPos;
			if ((size_t)cPos < src.size())
				dest[i] = src[(size_t)cPos] * (oldPos - fPos) + src[(size_t)fPos] * (cPos - oldPos);
			else
				dest[i] = src.back();
		}
	}
}


/// calculates covariance between vectors a and b
template<class T>
T covariance(const std::vector<T>& a, T aMean, size_t aStart, const std::vector<T>& b, T bMean, size_t bStart, size_t len) {
	T res(0);
	
	for (size_t i = aStart, j = bStart; i < aStart + len; ++i, ++j) {
		res += (a[i] - aMean)*(b[j] - bMean);
	}
	
	return res / T(len);
}


/// container of result for the convolution function
/// @warning the name should be changed since now it is used in correlation and other functions too
template<class T = double>
struct ConvolutionResult {
	T bestValue;
	int bestOffset;
	
	ConvolutionResult() : bestValue(0), bestOffset(0) {}
};


/// cross correlation between vectors a and b
/// if b is shorter than a, different offsets for b are tried, and only the one giving the best result is returned
template<class T> 
ConvolutionResult<T> crossCorrelation(const std::vector<T>& a, const std::vector<T>& b) {
	//std::cout << "cross-correlation of vector sizes: " << a.size() << ", " << b.size() << "\n";
	ConvolutionResult<T> res;
	
	using std::min;
	using std::max;
	
	// for (int offs = 1-b.size(); offs < (int)a.size()-1; ++offs) {
	for (int offs = 0; offs < int(a.size()-b.size()); ++offs) {
		double cor = 0;
		double sumSqrA = 0;
		for (int ia = max(0, offs); ia < (int)min(a.size(), b.size() + offs); ++ia) {
			cor += a[ia] * b[ia - offs];
			sumSqrA += a[ia] * a[ia];
		}
		
		cor /= sumSqrA;
		
		if ((cor > res.bestValue) || (offs == 1-(int)b.size())) {
			res.bestValue = cor;
			res.bestOffset = offs;
		}
	}
	
	return res;
}


/// Parsons coefficient of correlation between vectors a and b
template<class T> 
ConvolutionResult<T> statisticalCorrelationCoeff(const std::vector<T>& a, const std::vector<T>& b, int offs = -1) {
	// find the best offset
	ConvolutionResult<T> res;
	if (offs < 0) {
		res = crossCorrelation(a, b);
	} else {
		res.bestOffset = offs;
	}
	
	// calculate mean and std. dev. of the overlapping parts of both functions
	T meanA, stdDevA, meanB, stdDevB, covar;
	meanAndVar(meanA, stdDevA, a, std::max(0, res.bestOffset), std::min(a.size(), b.size() + res.bestOffset));
	stdDevA = sqrt(stdDevA);
	meanAndVar(meanB, stdDevB, b, std::max(0, -res.bestOffset), std::min(b.size(), a.size() - res.bestOffset));
	stdDevB = sqrt(stdDevB);
	covar = covariance(
		a, meanA, std::max(0, res.bestOffset), 
		b, meanB, std::max(0, -res.bestOffset), 
		std::min(a.size(), b.size() + res.bestOffset) - std::max(0, res.bestOffset));
	
	res.bestValue = covar / (stdDevA * stdDevB);
	
	/*
	std::cout << "covar = " << covar << "\n" << "stdDev = " << stdDevA << ", " << stdDevB << "\n" <<
		"mean = " << meanA << ", " << meanB << "\n" << "k = " << covar << "/" <<
		stdDevA * stdDevB << " = " << (covar / (stdDevA * stdDevB)) << "\n" << 
		"best offset is " << res.bestOffset << "\n";
	*/
	
	return res;
}


/// RMS (root mean squared) of the difference between vectors a and b
/// @param [in] offset is offset of b compared to a; if it is -1, then best offset is searched for using convolution
template<class T> 
ConvolutionResult<T> rms(const std::vector<T>& a, const std::vector<T>& b, int offs = -1) {
	// find the best offset
	ConvolutionResult<T> res;
	if (offs < 0) {
		res = crossCorrelation(a, b);
	} else {
		res.bestOffset = offs;
	}
	
	size_t len = std::min(a.size(), b.size() + res.bestOffset) - std::max(0, res.bestOffset);	
	std::vector<T> aMinusB(len);
	for (size_t i = 0; i < len; ++i)
		aMinusB[i] = a[std::max(0, res.bestOffset) + i] - b[std::max(0, -res.bestOffset) + i];
	aMinusB = sqr(aMinusB);
	T mean, stdDev;
	meanAndVar(mean, stdDev, aMinusB, 0, len);
	
	res.bestValue = sqrt(mean);
	
	return res;
}


/// correlation coefficient between vectors (fast to calculate, a bit different than statistical correlation)
/// @warning vector correlation is sensitive to the vectors' mean
template<class T> 
ConvolutionResult<T> vectorCorrelationCoeff(const std::vector<T>& a, const std::vector<T>& b, int offs = -1) {
	// find the best offset
	ConvolutionResult<T> res;
	if (offs < 0) {
		res = crossCorrelation(a, b);
	} else {
		res.bestOffset = offs;
	}

	// r = dot product of vectors divided by product of vector lengths
	size_t len = std::min(a.size(), b.size() + res.bestOffset) - std::max(0, res.bestOffset);	
	res.bestValue = 
		dot(a, std::max(0, res.bestOffset), b, std::max(0, -res.bestOffset), len) / 
		( sqrt(dot(a, std::max(0, res.bestOffset), a, std::max(0, res.bestOffset), len)) *
		  sqrt(dot(b, std::max(0, -res.bestOffset), b, std::max(0, -res.bestOffset), len)) );
		
	return res;
}


/// a function of similarity between two curves (represented as vectors of samples)
/// was used when U wave was sought, because if there are multiple peaks of variuos sizes, this 
/// function tries to fit them all with the same weight 
template<class T> 
ConvolutionResult<T> devFromLinear(const std::vector<T>& a, const std::vector<T>& b, T bOfs, int offs = -1) {
	static const double reallyBadFitness = 1e-30;
	// find the best offset
	ConvolutionResult<T> res;
	if (offs < 0) {
		res = crossCorrelation(a, b);
	} else {
		res.bestOffset = offs;
	}
	
	// r = dot product of vectors divided by product of vector lengths
	size_t len = std::min(a.size(), b.size() + res.bestOffset) - std::max(0, res.bestOffset);	
	std::vector<double> d(len);
	
	double aMin = 0, aMax = 0;
	minAndMax(aMin, aMax, a);
	double aMult = 1.0 / (aMax - aMin);
	double aOfs = bOfs; //1 - (aMin * aMult);
	
	for (size_t i = std::max(0, res.bestOffset), j = std::max(0, -res.bestOffset), k = 0; k < len; ++k, ++i, ++j)
		d[k] = (a[i] * aMult + aOfs) / (b[j] + bOfs);
	
	double mean, var;
	meanAndVar(mean, var, d, 0, len);
	
	if (var == 0)
		res.bestValue = 1/reallyBadFitness;
	else
		res.bestValue = sqrt(var); //mean / sqrt(var);
	
	return res;
}


#endif // VECTORMATH_H_INCLUDED
