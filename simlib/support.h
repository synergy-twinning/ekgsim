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

#ifndef SUPPORT_H_INCLUDED
#define SUPPORT_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <queue>
#include <map>
#include <cmath>
#include <ctime>

#include <Ini.h>
#include <Pattern/Hypermatrix.h>
#include "columnFile.h"
#include "Wohlfart.h"



	// ***********************************************************************************************
	// function filenameExtension
	//
	/// retrieves the part of file name after the last dot
	///
	template<class T>
	std::basic_string<T> filenameExtension(const std::basic_string<T>& filename) {
		for (int i = (int)filename.size()-1; i >= 0; --i) {
			if (filename[i] == '.') {
				return std::basic_string<T>(filename.begin() + i, filename.end());
			}
		}
		return "";
	}


namespace SimLib {

	using std::string;
	using std::cout;
	using std::cerr;
	using std::endl;
	using std::ofstream;
	using std::ifstream;
	using std::min;
	using std::max;
	using std::priority_queue;
	using std::getline;
	using std::count;
	using std::istringstream;
	using std::ostringstream;
	using std::swap;

	using Pattern::Region;


	// ***********************************************************************************************
	// support functions and classes
	// ***********************************************************************************************


	// ***********************************************************************************************
	// struct StreamRedirector
	//
	/// Used to redirect one of the default streams (cout, cerr, clog) into a file or another stream
	struct StreamRedirector {
		std::ostringstream messageBuffer;
		std::streambuf* oldBuffer;
		std::ostream& stream;

		/// supply the stream to be redirected to the ctor
		StreamRedirector(std::ostream& s) : stream(s) {
			oldBuffer = stream.rdbuf();
			stream.rdbuf(messageBuffer.rdbuf());
		}

		/// when going out of scope, redirect the stream back to its original destination
		~StreamRedirector() {
			stream.rdbuf(oldBuffer);
		}
	};


	#ifdef WIN32
		static const char dirSeparator = '\\';
	#else
		static const char dirSeparator = '/';
	#endif


	// ***********************************************************************************************
	// class ScopeTimer
	//
	/// ScopeTimer is used to time the time of execution of a piece of code in a single C++ scope (
    /// code that is between curly braces {}). This is only base class and only provides building blocks
    /// for building scope timers with specific actions in their dtors
	class ScopeTimer {
		clock_t								startTime;

	public:
		ScopeTimer() {
			startTime = clock();
		}

		void reset() {
			startTime = clock();
		}

		clock_t check() const {
			return clock() - startTime;
		}

		float secondsPassed() const {
			return (float)check() / (float)CLOCKS_PER_SEC;
		}
	};


	// ***********************************************************************************************
	// class ScopeTimerForLogging
	//
	/// Scope timer that logs the time passed when it reaches its dtor.
	class ScopeTimerForLogging : public ScopeTimer {
		std::ostream&						log;
	public:
        /// ctor requires log stream and message. Message is written to log stream immediatelly,
        /// time passed is then added to log stream when the timer goes out of scope
		ScopeTimerForLogging(std::ostream& logStream, const std::string& message) : log(logStream) {
			log << message;
		}

		~ScopeTimerForLogging() {
			log << "done (" << secondsPassed() << "s)\n";
		}

		/// reload can be used to force timer to output current time elapsed and to start a new timer
		/// together with a new log entry
		void reload(const std::string& message) {
			log << "done (" << secondsPassed() << "s)\n";
			reset();
			log << message;
		}
	};


	// ***********************************************************************************************
	// operator <<
	//
	/// output operator for Vector template class
	///
	template<class charT, class charTraits, class T, size_t N>
	std::basic_ostream<charT, charTraits>& operator<< (std::basic_ostream<charT, charTraits>& out, const Pattern::Vector<T, N>& v) {
		if (N > 0) {
			for (size_t i = 0; i < N-1; ++i)
				out << v[i] << ", ";
		}
		out << v[N-1];
		return out;
	}


	// ***********************************************************************************************
	// function sqr
	//
	template<class T>
	T sqr(const T var) {
		return (var * var);
	}


	// ***********************************************************************************************
	// function powOf
	//
	/// Returns parameter on the power that is specified as template parameter (unsigned integer).
	///
	template<size_t N>
	struct PowOf {
		template<class T>
		inline T operator() (const T var) {
			return var*PowOf<N-1>(var);
		}
	};


	template<>
	struct PowOf<1u> {
		template<class T>
		inline T operator() (const T var) {
			return var;
		}
	};


    /// interpolates two numbers given two weights between 0 and 1 (they don't need to add up to 1)
	template<class T>
	double interpolate(T a, T wa, T b, T wb) {
		return (a * wa + b * wb) / (wa + wb);
	}


	// ***********************************************************************************************
	// function apd90
	//
	/// APD90 calculates the 90% repolarization time, givem a vector of data
	/// @warning Wohlfart and WohlfartPlus have a better implementation of APD90 function integrated
	template<class Vec>
	double apd90(const Vec& v) {
		double target = *std::max_element(v.begin(), v.end()) * 0.1 + *std::min_element(v.begin(), v.end()) * 0.9;
		double time = -1.0;
		for (size_t i = 1; i < v.size(); ++i) {
			if ((v[i-1] > target) && (v[i] < target)) {
				time = interpolate(double(i-1), (v[i-1] - target),
					double(i), (target - v[i]));
			}
		}
		return time;
	}

} // namespace SimLib

#endif // SUPPORT_H_INCLUDED
