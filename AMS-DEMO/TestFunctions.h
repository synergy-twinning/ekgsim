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

#ifndef TESTFUNCTIONS_H_INCLUDED
#define TESTFUNCTIONS_H_INCLUDED


#include "Array.h"
#include <vector>
#include <cassert>


namespace TestFunctions {
	
	/// minimization test functions
	
	/// #################################################################################
	/// single dimensional (return value is a scalar)
	/// #################################################################################
	template<size_t D = 2>
	struct Sphere {
		typedef double Value;
		typedef Array<double, D> Input;
		typedef void Properties;
		enum {geneBounds = false};
		
		double operator() (const Input& solution, Value& value, double& violation, ...) {
			value = 0;
			for (size_t i = 0; i < D; ++i)
				value += (solution[i]+5) * (solution[i]+5);
			violation = 0.0;
			return violation;
		}
		
		inline void normalize(Input&) const {}
	};
	
	
	/// #################################################################################
	/// multidimensional (return value is a vector)
	/// #################################################################################
	
	/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	/// ZDT1
	///
	/// dimensionality D must be > 0
	///
	template<size_t D>
	struct ZDT1 {
		typedef Array<double, 2> Value;
		typedef Array<double, D> Input;
		typedef void Properties;
		enum {geneBounds = true};
		
		static const size_t dimensionality = D;
		
		void operator() (const Input& solution, Value& result, double& violation, ...) const {
			result[0] = solution[0];
			double sum = 1;
			for (size_t i = 1; i < solution.size(); ++i) 
				sum += 9 * solution[i] / (solution.size() - 1);
			result[1] = sum * (1 - sqrt(solution[0] / sum));
			
			violation = 0.0;
		}
		
		void normalize(Input& solution) const {
			for (size_t i = 0; i < D; ++i)
				solution[i] = std::min(std::max(solution[i], 0.0), 1.0);
		}
	};
	
}


#endif // TESTFUNCTIONS_H_INCLUDED
