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

#ifndef PARALLELFRAMEWORK_H_INCLUDED
#define PARALLELFRAMEWORK_H_INCLUDED


// this file should be used instead of jobDistributer



#include "mpiWrapper.h"
#include <iostream>
#include <algorithm>
#include <deque>
#include <map>
#include <memory>
#include <stdexcept>
#include "VectorArithmetics.h"
#include "timer.h"
#include "utilities.h"


// hidden namespace for message tags
namespace {
	struct Tag {
		enum {
			exit = 0,
			newJob,
			finishedJob,
			newJobAndInvalidateQueue,
			gatherer,
			invalid
		};
	};
}


/**
/// This is how Job class should look like (replacing int with any kind of type)
	struct JobTemplate {
		int operator() (int) {return 0;}
	};
**/


template<class Input>
class JobQueue {
public:
	struct JobInfo {
		Input input;
		int id;
		
		// lifeTimer times the total lifetype of a job (from its generation to the
		// return of its output) on the master side
		PrecisionTimer1 lifeTimer;
		
		JobInfo(const Input& in, int theId) : input(in), id(theId) {
			lifeTimer.start();
		}
	};
	
protected:
	Mpi::Communicator communicator;
	// new jobs are appended, while processed jobs are taken from deque front
	std::deque<JobInfo> jobs;
	// when allDone is true, worker can safely quit 
	// notice that empty jobs vector does not mean there will be no more jobs
	bool allDone;
	BinaryStream stream;
	
public:
	JobQueue(const Mpi::Communicator& comm) : communicator(comm), allDone(false) {}
	JobQueue(const JobQueue& other) : communicator(other.communicator), allDone(false) {}
	
	const JobQueue& operator=(const JobQueue& other) {communicator = other.communicator; return *this;}
	
	size_t size() const {
		return jobs.size();
	}
	
	bool finished() const {return (allDone && jobs.empty());}
	// master should call finish when there is no more jobs to do
	void finish() {allDone = true;}
	
	// add a job by specifying its input (called by master)
	void add(const Input& in, int id) {
		jobs.push_back(JobInfo(in, id));
	}
	
	// receive a job from a remote process (blocking function, so make sure
	// a job is waiting when calling it)
	void receiveRemoteJob(Mpi::Status& status) {
		Input in;
		int id;
		stream.clear();
		Mpi::receive(stream, status, communicator);
		stream >> in >> id;
		add(in, id);
//		std::cerr << Mpi::Communicator().getRank() << " received   " << id << ",   in queue: ";
//		for (size_t i = 0; i < jobs.size(); ++i)
//			std::cerr << jobs[i].second << " ";
//		std::cerr << std::endl;
	}
	
	// called by master, this function inserts job into local queue and also
	// sends it to remote queue copy
	// when the remote queue processes this job, it should signal it and process
	// should manually remove it from local (this) queue by calling processNext
	void addRemote(const Input& in, int dest, int id) {
		static Mpi::Request request;
		stream.clear();
		stream << in << id;
		Mpi::send(stream, dest, Tag::newJob, request, communicator);
		add(in, id);
		request.wait();
	}
	
	JobInfo processNext() {
		JobInfo temp = jobs.front();
//		std::cerr << Mpi::Communicator().getRank() << " processing " << temp.second << ", from queue: ";
//		for (size_t i = 0; i < jobs.size(); ++i)
//			std::cerr << jobs[i].second << " ";
//		std::cerr << std::endl;
		jobs.pop_front();
//		std::cerr << Mpi::Communicator().getRank() << " removing job from local queue, size = " << size() << "\n";
		return temp;
	}
	
	// function called by workers, returns true if queue is not empty
	bool workerCheck(bool block = false) {
		int dummy;
		
		Mpi::Status status;
		if (block)
			status.wait(communicator);
		// process all waiting messages
		while (status.probe(communicator)) {
			switch (status.tag()) {
			case Tag::exit:
				Mpi::receive(dummy, status, communicator);
				finish();
				break;
			case Tag::newJob:
				receiveRemoteJob(status);
				break;
			default:
				throw std::runtime_error("unexpected tag received by worker!");
			}
		}
		return !jobs.empty();
	}
};


namespace {
	
	struct DummyControlJob {
	};
	
}


template<class Input, class Output, class Queue = JobQueue<Input> >
class PFramework {
	const Mpi::Communicator communicator;
	Output lastOut;
	bool purgeQueues; // not yet implemented
	// when blockingWorkerWait is true, workers will block inside the execute function, 
	// if the queue is empty (waiting for job updates), otherwise workers will simply
	// return true on empty queue
	bool blockingWorkerWait;
	
public:
	Queue queue;
	std::vector<Queue> remoteQueues;
	int targetQueueSize; // min size, try to always have at least this full
	int maxQueueSize;	 // not yet fully implemented (the algorithm does not require this value)
	BinaryStream stream;
	int masterRank;
	bool masterIsWorker;
	bool debugOutput;
	
public:
	PFramework(const Mpi::Communicator& com) : 
		communicator(com), 
		purgeQueues(false),
		blockingWorkerWait(true),
		queue(com),
		targetQueueSize(2),
		maxQueueSize(4),
		masterRank(0),
		masterIsWorker(false),
		debugOutput(true)
	{
		if (master()) {
			remoteQueues.resize(com.getSize(), queue);
			if (com.getSize() == 1)
				masterIsWorker = true;
		}
	}
	
	void setBlockingWorker(bool block = true) {
		blockingWorkerWait = block;
	}
	
	// create new id (every job gets its own id)
	int newId() {
		static int lastId = 0;
		return ++lastId;
	}
	
	// safest not to call this after the queues are in use (after execute has been 
	// called at least once (although it might at worst cause some workers to wait idle)
	void setQueueSizes(int target, int maximum) {
		targetQueueSize = target;
		maxQueueSize = maximum;
	}
	
	// fuction that differentiates between master and worker processes
	// returns false when execution is finished (if return is true excute should be
	// called again)
	// if worker blocking is enabled, workers will block inside the function, waiting 
	// for either the queue to refill or to receive exit signal.
	template<class Job, class ControlJob>
	bool execute(Job& job, ControlJob& controlJob, PrecisionTimer* communicationTimer = 0, PrecisionTimer* idleTimer = 0) {
		int rank = communicator.getRank();
		
		if (master()) {
			// synchronize local queue copies with remote queues
			synchronizeQueues(controlJob, communicationTimer);
			
			// count number of jobs in queues, sort queues by their size
			size_t totalQueuedJobs = 0;
			std::vector<size_t> queueSizes(remoteQueues.size(), 1);
			std::vector<size_t> queueIndices;
			
			if (purgeQueues) {
				queueIndices.resize(queueSizes.size());
			} else {
				for (size_t i = 0; i < remoteQueues.size(); ++i) {
					// if master is not worker, treat its queue as full
					if (((int)i == masterRank) && !masterIsWorker)
						queueSizes[i] = maxQueueSize;
					else
						queueSizes[i] = remoteQueues[i].size();
					totalQueuedJobs += queueSizes[i];
				}
				getSortedIndices(queueSizes, queueIndices);
			}
			
			// add jobs to queues, to keep them between the min and max sizes
			// (always target max size and equal sizes)
			// first determine number of desired jobs
			int numOfJobs = maxQueueSize*queueSizes.size() - totalQueuedJobs;
			// then let controller decide how many will be created
			numOfJobs = controlJob.moreJobs(numOfJobs);
			if (debugOutput)
                std::cerr << queueSizes << " + " << numOfJobs << " jobs -> ";
			while (numOfJobs > 0) {
				if (queueSizes.size() > 1) {
					// find first que (from sorted queues) that is shorter then the one before 
					//    it; queues were sorted from shortest to longest, therefore in first
					//    iteration, first queue will be found
					size_t shortestQueueIndex = 0;
					for (; shortestQueueIndex < (queueSizes.size() - 1); ++shortestQueueIndex) {
						size_t j = queueIndices[shortestQueueIndex];
						size_t j1 = queueIndices[shortestQueueIndex+1];
						if (queueSizes[j] < queueSizes[j1])
							break;
					}
					shortestQueueIndex = queueIndices[shortestQueueIndex];
					
					// create new job and add it to the selected queue;
					// distinguish between local queue and remote queues
					int id = newId();
					++queueSizes[shortestQueueIndex];
					if (rank != (int)shortestQueueIndex) {
						ScopeTimer t(*communicationTimer);
						remoteQueues[shortestQueueIndex].addRemote(controlJob.generateJob(id), shortestQueueIndex, id);
					} else
						remoteQueues[shortestQueueIndex].add(controlJob.generateJob(id), id);
				} else {
					int id = newId();
					remoteQueues[0].add(controlJob.generateJob(id), id);
				}
				--numOfJobs;
				++totalQueuedJobs;
			}
			if (debugOutput)
                std::cerr << queueSizes << std::endl;
			
			if (!masterIsWorker) {
				// repair totalQueuedJobs variable
				totalQueuedJobs -= queueSizes[rank];
			}
			
			// work on a job locally
			if(masterIsWorker && (remoteQueues[rank].size() > 0)) {
				typename Queue::JobInfo out = remoteQueues[rank].processNext();
				PrecisionTimer evaluationTimer;
				{
					ScopeTimer tt(evaluationTimer);
					lastOut = job(out.input);
				}
				out.lifeTimer.pause();
				controlJob.newResult(lastOut, out.id, rank, evaluationTimer.totalSeconds(), out.lifeTimer.totalSeconds());
				--totalQueuedJobs;
			}
			
			// controller says there will be absolutely no more jobs and there are no 
			// pending jobs -> tell everyone to wrap up
			if ((numOfJobs < 0) && (totalQueuedJobs == 0)) {
				finishAll(communicationTimer);
				return false;
			}
			
			// blocking wait if there is no more work to be done (and the blocking 
			// wait is enabled) and some workers are still busy
			if ((!masterIsWorker || (numOfJobs <= 0)) && blockingWorkerWait && 
				(remoteQueues[rank].size() == 0) && (totalQueuedJobs > 0)) 
			{
				blockingWait(idleTimer);
			}
			
			return true;
		} else { // worker
			ScopeTimer t(*communicationTimer);
			
			// early bail
			if (queue.finished())
				return false;
			
			// determine if there is a job available
			// also wait for job updates if blocking workers are enabled 
			bool jobAvailable = queue.workerCheck();
			if (blockingWorkerWait && !queue.finished() && !jobAvailable) {
				ExcludeScopeFromTimer t1(*communicationTimer);
				ScopeTimer t2(*idleTimer);
				if (debugOutput)
                    std::cerr << "[" << rank << "] is waiting\n";
				jobAvailable = queue.workerCheck(true);
			}
			/*
			// old blocking code
			while (blockingWorkerWait && !queue.finished() && !jobAvailable) {
				if (queue.workerCheck(communicator, true)) {
					jobAvailable = true;
					break;
				}
			}
			*/
			
			// work on a job and signal back results
			if (jobAvailable) {
				typename Queue::JobInfo nextIn = queue.processNext();
				PrecisionTimer evaluationTimer;
				{
					ScopeTimer tt(evaluationTimer);
					ExcludeScopeFromTimer et(*communicationTimer);
					lastOut = job(nextIn.input);
				}
				processResult(lastOut, nextIn.id, evaluationTimer.totalSeconds());
				
//				std::cerr << "sent " << nextIn.second << ";  " << lastOut.chromosome << " | " << lastOut.criteria << std::endl;
			}
			return !queue.finished();
		}
	}
	
	void processResult(const Output& result, int id, double evalTime) {
		stream.clear();
		stream << result << id << evalTime;
		Mpi::Request request;
		Mpi::send(stream, masterRank, Tag::finishedJob, request);
		request.wait();
	}
	
	// return true if process is master (rank == 0)
	bool master() const {
		return communicator.getRank() == masterRank;
	}
	
	// return last output (output is only defined after the first job processes)
	const Output& lastOutput() const {
		return lastOut;
	}
	
	// called by master process, function signals all processes to finish
	// it is a blocking command - blocks until every process accepts finish command
	void finishAll(PrecisionTimer* communicationTimer) {
		assert(master());
		
		ScopeTimer t(*communicationTimer);
		queue.finish();
		int size = communicator.getSize();
		std::vector<Mpi::Request> requests(size);
		
		for (size_t i = 1; (int)i < size; ++i) {
			int dummy;
			Mpi::send(dummy, i, Tag::exit, requests[i]);
		}
		
		for (size_t i = 1; (int)i < size; ++i) {
			requests[i].wait();
		};
	}
	
	// function called by master when no jobs are to be processed locally and some jobs
	// are in remote queues
	void blockingWait(PrecisionTimer* waitTimer = 0) {
		ScopeTimer t(*waitTimer);
		Mpi::Status status;
		
		// MPI_probe is an active awit - we do not want that, therefore
		// go to short sleep while there is no messages (just yield cpu to other processes)
//		PrecisionTimer tt; tt.start();
		size_t loops = 0;               // just for debugging active loop
		for (; !status.probe(communicator); ++loops) {
            timespec requestTime;
            requestTime.tv_nsec = 1;
            requestTime.tv_sec = 0;
            nanosleep(&requestTime, &requestTime);
		}
//        tt.pause(); std::cout << "waited " << tt.totalSeconds() << "s; made " << loops << "loops\n";
//		status.wait(communicator);
	}
	
	// function called by master, returns true if queue is not finished
	template<class ControlJob>
	void synchronizeQueues(ControlJob& controlJob, PrecisionTimer* communicationTimer) {
		assert(master());
		ScopeTimer t(*communicationTimer);
		
		Mpi::Status status;
		// process all waiting messages
		while (status.probe(communicator)) {
			stream.clear();
			Mpi::receive(stream, status);
			if (status.tag() == Tag::finishedJob) {
				int id;
				double evalTime;
				stream >> lastOut;
				stream >> id;
				stream >> evalTime;
//				std::cerr << "picked " << id << ";  " << lastOut.chromosome << " | " << lastOut.criteria << std::endl;
				// remove job from local queue
				{
                    ExcludeScopeFromTimer et(*communicationTimer);
                    typename Queue::JobInfo j = remoteQueues[status.source()].processNext();
                    j.lifeTimer.pause();
                    controlJob.newResult(lastOut, id, status.source(), evalTime, j.lifeTimer.totalSeconds());
                    if ((id != j.id) && debugOutput) {
                        std::cerr << "warning: in synchronization, individual " << j.id
                            << " from remote queue replaced individual " << id
                            << " from local queue " << status.source() << "\n"; 
                    }
				}
			} else {
			}
		}
	}
};


class Gatherer {
	const Mpi::Communicator comm;
	int countdownReceived;
	
public:
	Gatherer(const Mpi::Communicator& c = Mpi::Communicator()) : comm(c), countdownReceived(c.getSize()-1) {}
	
	// returns true when gathering is finished for the process that called this function
	// (for master it is finished as soon as it receives data from all the processes and
	// for slaves it is over as soon as they send the data)
	template<class Receiver, class Sender>
	bool gatherIn(Receiver& receiver, const Sender& sender, int destinationRank) {
		if (comm.getRank() == destinationRank) {
			if (countdownReceived == 0)
				return true;
			// gatherer (master)
			Mpi::Status status;
			BinaryStream stream;
			while (status.probe(comm)) {
				switch (status.tag()) {
				case Tag::gatherer:
					Mpi::receive(stream, status, comm);
					--countdownReceived;
					if (!receiver.receive(status.source(), stream) || (countdownReceived < 0))
						throw std::runtime_error("error in gatherer");
					stream.clear();
					break;
				default:
					throw std::runtime_error("detected an invalid message during gathering");
					break;
				}
			}
			return (countdownReceived == 0);
		} else {
			// sender (slave)
			Mpi::Request request; 
			Mpi::send(sender.sendBuffer(), destinationRank, Tag::gatherer, request, comm);
			request.wait();
			return true;
		}
	}
};

#endif // PARALLELFRAMEWORK_H_INCLUDED
