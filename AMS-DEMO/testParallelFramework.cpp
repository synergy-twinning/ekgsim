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
