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

#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED


#include <vector>
#include <cstdio> // for size_t 
#include <ctime>
#include <cerrno>
#include <stdexcept>


#ifdef WIN32
	#include <windows.h>
	#include <time.h>
	#define time_type LARGE_INTEGER
	#define MAC_getTime(X) QueryPerformanceCounter(&X);

	#define MAC_addDifference(dest, a, b) dest.QuadPart += a.QuadPart - b.QuadPart;
#else
	#include <sys/time.h>
	#define time_type timeval
	#define MAC_getTime(X) gettimeofday(&X,0)

	#define MAC_addDifference(dest, a, b) dest.tv_sec += a.tv_sec - b.tv_sec; dest.tv_usec += a.tv_usec - b.tv_usec;
#endif //WIN32


namespace {

	double toSeconds(const time_type& totalTime) {
#ifdef WIN32
		LARGE_INTEGER freq;
		QueryPerformanceFrequency(&freq);
		return (double)totalTime.QuadPart / (double)freq.QuadPart;
#else
		return (double)totalTime.tv_sec + totalTime.tv_usec*0.000001;
#endif //WIN32
	}

	double toMicroSeconds(const time_type& totalTime) {
#ifdef WIN32
		LARGE_INTEGER freq;
		QueryPerformanceFrequency(&freq);
		return (double)totalTime.QuadPart * 1000000.0 / (double)freq.QuadPart;
#else
		return (double)totalTime.tv_sec*1000000.0 + totalTime.tv_usec;
#endif //WIN32
	}
#ifdef WIN32
    void interruptableSleep(size_t nanos) {
        Sleep(0);
	}
#else
	void interruptableSleep(size_t nanos) {
	    timespec requestTime;
		requestTime.tv_nsec = nanos % 1000000;
		requestTime.tv_sec = nanos / 1000000;
        while(nanosleep(&requestTime, &requestTime)==-1 && errno == EINTR);
	}
#endif //WIN

}


class PrecisionTimer1 {
	time_type startTime;
	time_type totalTime;

public:
	PrecisionTimer1() {
		time_type temp = {0};
		totalTime = temp;
	}

	void start() {
		MAC_getTime(startTime);
	}

	time_type read() {
		time_type temp;
		MAC_getTime(temp);
		return temp;
	}

	void pause() {
		time_type endTime;
		MAC_getTime(endTime);
		MAC_addDifference(totalTime, endTime, startTime);
	}

	double totalSeconds() const {
		return toSeconds(totalTime);
	}

	operator time_type() {
		return totalTime;
	}
};


class PrecisionTimerSequence {
	std::vector<std::pair<time_type, int> > sequence;

public:
	void reserve(size_t rsize) {
		sequence.reserve(rsize);
	}

	void addTime(const time_type& t, int id) {
		sequence.push_back(std::make_pair(t, id));
	}

	size_t size() const {
		return sequence.size();
	}

	std::pair<double, int> at(size_t i) const {
		time_type temp = {0};
		MAC_addDifference(temp, sequence[i].first, sequence[0].first);
		return std::make_pair(toSeconds(temp), sequence[i].second);
	}
};


class PrecisionTimer {
	time_type startTime;
	time_type totalTime;
	PrecisionTimerSequence* seq;
	int id;
	bool paused;

public:
	PrecisionTimer() : seq(0), paused(true) {
		time_type temp = {0};
		totalTime = temp;
	}

	PrecisionTimer(PrecisionTimerSequence& s, int i) : seq(&s), id(i), paused(true) {
		time_type temp = {0};
		totalTime = temp;
	}

	void start() {
	    if (paused) {
            paused = false;
            MAC_getTime(startTime);
            if (seq)
                seq->addTime(startTime, id);
	    } else
            throw std::runtime_error("not paused when trying to start");
	}

	time_type read() {
		time_type temp;
		MAC_getTime(temp);
		return temp;
	}

	void pause() {
	    if (!paused) {
	        paused = true;
            time_type endTime;
            MAC_getTime(endTime);
            MAC_addDifference(totalTime, endTime, startTime);
            if (seq)
                seq->addTime(endTime, -id);
	    } else throw std::runtime_error("already paused when trying to pause");
	}

	double totalSeconds() const {
	    if (!paused)
            throw std::runtime_error("timer not paused when total number of seconds is being read!");
		return toSeconds(totalTime);
	}

	operator time_type() {
		return totalTime;
	}

	void clear() {
	    paused = true;
		time_type temp = {0};
		totalTime = temp;
	}

	/// wait until total time passed is larger than or equal to "tm"
	void waitTotal(double tm) {
	    //do {MAC_getTime(t); MAC_addDifference(tt, t, startTime);} while (toMicroSeconds(tt) < tm);
	    do {
	        time_type t, tt = {0};
	        MAC_getTime(t);
	        MAC_addDifference(tt, t, startTime);
	        if (toSeconds(tt) >= tm*0.000001)
                break;
	    } while(true);
	}
};


// create a scope timer to add the duration of current scope to selected PrecisionTimer
class ScopeTimer {
	PrecisionTimer& pt;

public:
	ScopeTimer(PrecisionTimer& t) : pt(t) {
		if (&pt != 0) pt.start();
	}

	~ScopeTimer() {
		if (&pt != 0) pt.pause();
	}
};


class ExcludeScopeFromTimer {
	PrecisionTimer& pt;

public:
	ExcludeScopeFromTimer(PrecisionTimer& t) : pt(t) {
		if (&pt != 0) pt.pause();
	}

	~ExcludeScopeFromTimer() {
		if (&pt != 0) pt.start();
	}
};



#endif // TIMER_H_INCLUDED
