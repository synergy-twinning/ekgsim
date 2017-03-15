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

#ifndef NONLINEAR_FIT_H_INCLUDED
#define NONLINEAR_FIT_H_INCLUDED

/** @file nonlinearFit.h
    @brief optimization function for curve fitting
**/

#include "vectorMath.h"
// typedef std::vector<double> TFunction;


/** @fn template<class Func, class Vec> int steepestDescend(Func& f, Vec& x0, const Vec& d, double stepSize = 0.5, double epsilon = 1e-8, int iterations = 100) {
    optimization function based on the steepest descend principle

    @param [in] stepSize defines the starting step size, which goes down by factor 0.5, when the algorithm  discovers that it is too large for further iterations
    @param [in] epsilon algorithm searches until the change in solutions is smaller than epsilon
    @param [in] iterations specifies the maximum number of iterations
    @param [in] f price function used in optimization
    @param [in,out] x0 initial solution
    @param [in] d offset used to calculate the gradient of f (grad(f, x) = (f(x+d) - f(x)) / d)

    Starting with a guess x0, optimize function towards target value through
    the method of the steepest descend (gradient is calculated numerically
    through an offset d from a given point). Result is stored in x0

    @todo test optimization in a controlled environment (not on APs)
**/
template<class Func, class Vec>
int steepestDescend(Func& f, Vec& x0, const Vec& d, double stepSize = 0.5, double epsilon = 1e-8, int iterations = 100) {
	// calculate gradient in x0
	double dLen = 0.0;
	for (size_t i = 0; i < d.size(); ++i) {
		dLen += d[i] * d[i];
	}
	dLen = sqrt(dLen);

	Vec grad = x0;
	Vec oldGrad = grad;
	for (size_t i = 0; i < oldGrad.size(); ++i)
		oldGrad[i] = 0;

	double y0 = f(x0);
//	std::cout.precision(2);
//	std::cout.width(7);
//	std::cout.fill(' ');
//	std::cout.setf(std::ios::scientific, std::ios::floatfield);
//	std::cout << "y0 = " << y0 << "\n";

	Vec move = d;

	for (; (iterations > 0) && (y0 > epsilon); --iterations) {
//		std::cout << std::scientific;
		// calculate gradient in x0
//		std::cout << "grad: ";
		double gradLen = 0;
		bool stepChange = false;
		for (size_t i = 0; i < x0.size(); ++i) {
			Vec x1 = x0;
			if (d[i] != 0) {
				x1[i] += d[i] * .001;
				grad[i] = (f(x1) - y0) / (d[i] * .001);
				gradLen += grad[i] * grad[i];
				if (grad[i] * oldGrad[i] < 0) {
					// change in gradient sign, step is too large
					move[i] *= 0.5;
					stepChange = true;
				} else if (fabs(grad[i]) > 0.75*fabs(oldGrad[i])) {
					move[i] *= 1.5;
				}
//				std::cout << grad[i] << " ";
			} else {
				grad[i] = 0;
			}
		}
//		std::cout << "; \t";
		gradLen = stepSize / sqrt(gradLen);

		// is gradient ok? (is there no change of sign in gradient)
		oldGrad = grad;

		Vec x1 = x0;
		// move in the direction of gradient
		for (size_t i = 0; i < x0.size(); ++i) {
			x1[i] -= stepSize * ((grad[i] > 0) ? move[i] : -move[i]);
		}
		double y1 = f(x1);

		if (y1 < y0) {
			y0 = y1;
			x0 = x1;
		} else {
			if (!stepChange)
				stepSize *= 0.5;
		}

//		std::cout.precision(4);
//		std::cout << "y0 = " << std::fixed << y0 << ", stepSize=" << stepSize << "\n";
	}

//	std::cout.precision(-1);
//	std::cout.unsetf(std::ios::fixed);

	return iterations;
}


#endif // NONLINEAR_FIT_H_INCLUDED
