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

#include "ParallelFramework.h"
#include "pDemo.h"
#include "TestFunctions.h"

#include <iostream>
#include <algorithm>
#include <windows.h>


struct Complex {
	std::vector<int> vec;
	double d;
	int i;
};


BinaryStream& operator<< (BinaryStream& bs, const Complex& c) {
	return bs << c.vec << c.d << c.i;
}


BinaryStream& operator>> (BinaryStream& bs, Complex& c) {
	bs >> c.vec;
	bs >> c.d;
	bs >> c.i;
	return bs;
}


struct Job {
	Complex operator() (Complex& in) {
		std::cerr << "job " << in.i << " is being executed on cpu " << Mpi::Communicator().getRank() << "\n";
		Sleep(0); //forfeit cpu time
		in.d = -in.d;
		in.vec.resize(4);
		in.vec[0] = Mpi::Communicator().getRank();
		return in;
	}
};


struct Bureaucracy {
	int jobId;
	
	Bureaucracy() : jobId(0) {}
	
	Complex generateJob() {
		Complex j;
		std::cerr << "generating new job " << jobId << "\n";
		if (jobId > 19)
			throw "job id overflow";
		j.i = jobId++;
		j.vec.resize(4, 0);
		j.d = 2+j.i*0.031;
		return j;
	}
	
	// return the number of jobs that can be generated at the call time;
	// supply the function with an ideal number of jobs to be generated
	// return value <0 will be a signal that generation has stopped
	// return value 0 will signal that no job can be generated at this time but will be 
	// in the future
	int moreJobs(int preferredNum) {
		if (jobId == 20)
			return -1;
		else
			return std::min(20 - jobId, preferredNum);
	}
	
	void newResult(Complex res) {
		std::cerr << "recieved output " << res.i << "," << res.d << "," << res.vec[0] << "\n";
	}
};


int main(int argc, char** argv) {
	try {
		Mpi::Environment mpi(argc, argv);
		Mpi::Buffer buffer(10000);
		
		/*
		// test mpi I communication
		int from = (mpi.getCommunicator().getRank() + mpi.getCommunicator().getSize() - 1) % mpi.getCommunicator().getSize();
		int to = (mpi.getCommunicator().getRank() + 1) % mpi.getCommunicator().getSize();
		int tag = 1;
		
		if (mpi.getCommunicator().getRank() == 0) {
			int data;
			Mpi::Status status;
			for (int i = 1; i < mpi.getCommunicator().getSize(); ++i) {
				while (!status.probe())
					std::cout << "\r waiting " << i << "\r";
				Mpi::receive(data, status.source(), 0);
				std::cout << "\r recieved " << data << " from " << status.source() << "\n";
			}
		} else {
			Mpi::Request dummy;
			
			int data = mpi.getCommunicator().getRank();
			Mpi::send(data, 0, 0, dummy);
			std::cout << "\r sent " << data << "\n";
		}
		*/
		
		// test parallel framework
		size_t numProcesses = mpi.getCommunicator().getSize();
		Job job;
		Bureaucracy bureau;
		PFramework<Complex, Complex> fw(mpi.getCommunicator());
		
		if (fw.master()) {
			while (fw.execute(job, bureau));
			std::cerr << "master finished\n";
		} else {
			while (fw.execute(job, bureau));
			std::cerr << "worker " << Mpi::Communicator().getRank() << " finished\n";
		}
	} catch (...) {
		std::cout << "exception!\n";
	}
	return 0;
}
