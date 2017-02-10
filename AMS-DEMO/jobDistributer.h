#ifndef JOB_DISTRIBUTER_H_INCLUDED
#define JOB_DISTRIBUTER_H_INCLUDED


#include "mpiWrapper.h"
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "VectorArithmetics.h"
#include "timer.h"

/// This is how Job class should look like
/**
	struct Job {
		int operator() (int) {return 0;}
	};
}
**/

template<class T>
inline size_t size(T) {
	return 1;
}

template<class T>
inline size_t size(const std::vector<T>& v) {
	return v.size();
}


template<class T>
inline void resize(T, size_t) {}

template<class T>
inline void resize(std::vector<T>& v, size_t newSize) {
	v.resize(newSize);
}

template<class T, size_t N>
inline void resize(Array<T, N>& v, size_t newSize) {
	assert(N == newSize);
}


template<class T>
inline T* begin(T& var) {
	return &var;
}

template<class T>
inline typename std::vector<T>::iterator begin(std::vector<T>& vec) {
	return vec.begin();
}

template<class T, size_t N>
inline T* begin(Array<T, N>& vec) {
	return vec.begin();
}


template<class T>
inline T* end(T& var) {
	return begin(var)+1;
}

template<class T>
inline typename std::vector<T>::iterator end(std::vector<T>& vec) {
	return vec.end();
}

template<class T, size_t N>
inline T* end(Array<T, N>& vec) {
	return vec.end();
}


template<class Chromosome, class Criteria>
class JobDistributer {
	const Mpi::Communicator 	communicator;
	
public:
	JobDistributer(const Mpi::Communicator& com) : communicator(com) {}

	Mpi::Communicator comm() {return communicator;}
	
	template<class Job>
	void executeWorker(Job& job, PrecisionTimer* pt = 0) {
		ScopeTimer t(*pt);
		
		int n = 0;
		size_t inBufferSize;
		size_t outBufferSize;
		Mpi::receive(n, 0, 1);
		
		Chromosome in;
		for (int i = 0; i < n; ++i) {
			BinaryStream inStream;
			Mpi::receive(inStream, 0, 1);
			inStream >> in;
//			std::cout << "worker " << communicator.getRank() << " received " << (i+1) << " of " << n << std::endl;
			Criteria out;
			double violation;
			{
				ExcludeScopeFromTimer et(*pt);
				violation = job(in, out);
			}
//			std::cout << "worker " << communicator.getRank() << " executed job" << std::endl;
			
			BinaryStream outStream;
			outStream << in << out << violation;
			Mpi::send(outStream, 0, 1);
//			std::cout << "worker " << communicator.getRank() << " sent     " << (i+1) << " of " << n << " (" << outBufferSize << ")" << std::endl;
		}
	}
	
	template<class InputCont, class Job>
	void executeBoss(InputCont& in, Job& job, PrecisionTimer* pt = 0) {
		ScopeTimer t(*pt);
		
		int nCPU = communicator.getSize();
		int nThis = communicator.getRank();
		if (nThis != 0)
			throw std::runtime_error("distributer node can only have rank 0!");
		
		if (nCPU == 1) {
			// no parallelism
			ExcludeScopeFromTimer et(*pt) ;
			for (size_t i = 0; i < in.size(); ++i)
				in[i].violation = job(in[i].chromosome, in[i].criteria);
		} else {
			// count number of individuals that need evaluation (violation must be -1)
			size_t numEvals = 0;
			for (size_t i = 0; i < in.size(); ++i)
				if (in[i].violation == -1) ++numEvals;
			
			size_t sendSize = size(in[0].chromosome);
			size_t receiveSize = size(in[0].chromosome) + size(in[0].criteria) + 1;
			
			// send orders (num of jobs to execute, size of chromosome and size of 
			// returned buffer - the whole individual [chromosome vector, criteria 
			// vector and violation code])
			for (int n = 1; n < nCPU; ++n) {
//				std::cout << "boss sent order " << (in.size() + nCPU - 1 - n) / (nCPU - 1) << " to " << n << std::endl;
				Mpi::send((int)((numEvals + nCPU - 1 - n) / (nCPU - 1)), n, 1);
				Mpi::send(sendSize, n, 1);
				Mpi::send(receiveSize, n, 1);
			}
			
			std::vector<double> receiveBuffer(receiveSize);
			
			// send data
			// indices hold index of each evaluated individual for each CPU
			std::vector<size_t> indices;
			indices.reserve(nCPU - 1);
			for (size_t i = 0; i < in.size();) {
				indices.clear();
				for (int n = 1; (n < nCPU) && (i < in.size());) {
					if (in[i].violation == -1) {
						indices.push_back(i);
						++n;
					}
					++i;
				}
				
				for (int n = 1; (n < nCPU) && (n <= indices.size()); ++n) {
					BinaryStream inStream;
					inStream << in[indices[n-1]].chromosome;
					Mpi::send(inStream, n, 1);
//					std::cout << "boss sent     " << (i+n) << std::endl;
				}
				
				for (int n = 1; (n < nCPU) && (n <= indices.size()); ++n) {
					BinaryStream stream;
//					std::cout << "boss expecting to receive results from " << (i+n) << " (" << receiveSize << ") \n";
					Mpi::receive(stream, n, 1);
					stream >> in[indices[n-1]].chromosome >> in[indices[n-1]].criteria 
						>> in[indices[n-1]].violation;
//					std::cout << "boss received " << (i+n) << std::endl;
				}
			}
//			std::cout << "boss: sending finished" << std::endl;
		}
	}
	
	template<class InputCont, class Job>
	void executeBossAndWorker(InputCont& in, Job& job, PrecisionTimer* pt = 0) {
		ScopeTimer t(*pt);
		
		int nCPU = communicator.getSize();
		int nThis = communicator.getRank();
		if (nThis != 0)
			throw std::runtime_error("distributer node can only have rank 0!");
		
		if (nCPU == 1) {
			// no parallelism
			ExcludeScopeFromTimer et(*pt);
			for (size_t i = 0; i < in.size(); ++i)
				in[i].violation = job(in[i].chromosome, in[i].criteria);
		} else {
			// count number of individuals that need evaluation (violation must be -1)
			size_t numEvals = 0;
			for (size_t i = 0; i < in.size(); ++i)
				if (in[i].violation == -1) ++numEvals;
			
			// packet sizes
			size_t sendSize = in[0].chromosome.size();
			size_t receiveSize = in[0].chromosome.size() + size(in[0].criteria) + 1;
			
			// send orders (num of jobs to execute, size of chromosome and size of 
			// returned buffer - the whole individual [chromosome vector, criteria 
			// vector and violation code])
			for (int n = 1; n < nCPU; ++n) {
//				std::cout << "boss sent order " << (in.size() + nCPU - 1 - n) / (nCPU - 1) << " to " << n << std::endl;
				Mpi::send((int)((numEvals + nCPU - 1 - n) / nCPU), n, 1); //++++
				Mpi::send(sendSize, n, 1);
				Mpi::send(receiveSize, n, 1);
			}
			
			std::vector<double> receiveBuffer(receiveSize);
			
			// send data
			// indices hold index of each evaluated individual for each CPU
			std::vector<size_t> indices;
			indices.reserve(nCPU);
			for (size_t i = 0; i < in.size();) {
				// select only individuals with no violation
				indices.clear();
				for (int n = 0; (n < nCPU) && (i < in.size());) {
					if (in[i].violation == -1) {
						indices.push_back(i);
						++n;
					}
					++i;
				}
				
				// send individuals to evaluate
				for (int n = 1; (n < nCPU) && (n < indices.size()); ++n) {
					BinaryStream inStream;
					inStream << in[indices[n]].chromosome;
					Mpi::send(inStream, n, 1);
//					std::cout << "boss sent     " << (i+n) << std::endl;
				} 
				
				// execute one evaluation locally
				{
					ExcludeScopeFromTimer et(*pt);
					in[indices[0]].violation = job(in[indices[0]].chromosome, in[indices[0]].criteria);
				}
				
				// gather results
				for (int n = 1; (n < nCPU) && (n < indices.size()); ++n) {
					BinaryStream stream;
//					std::cout << "boss expecting to receive results from " << (i+n) << " (" << receiveSize << ") \n";
					Mpi::receive(stream, n, 1);
					stream >> in[indices[n-1]].chromosome >> in[indices[n-1]].criteria 
						>> in[indices[n-1]].violation;
				}
			}
//			std::cout << "boss: sending finished" << std::endl;
		}
	}
};


#endif // JOB_DISTRIBUTER_H_INCLUDED
