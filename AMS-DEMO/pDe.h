#ifndef DE_H_INCLUDED
#define DE_H_INCLUDED


#include "DeBase.h"
#include "ParallelNumericOptimizer.h"

#include <map>


/// ##################################################################################
/// template class PDe
/// used as an Algorithm policy for NumericOptimizer<...>;
///
template<class Chromosome, class Evaluation>
class PDe : 
	// inherit base DE definitions
	public DeBase<IndividualStruc<Chromosome, typename Evaluation::Value, typename Evaluation::Properties> >, 
	// make PDe conform the form of ParallelNumericOptimizer
	public ParallelNumericOptimizer<PDe<Chromosome, Evaluation>, Chromosome, Evaluation> 
{
public:
	friend class ParallelNumericOptimizer<PDe<Chromosome, Evaluation>, Chromosome, Evaluation>;
	typedef typename ParallelNumericOptimizer<PDe<Chromosome, Evaluation>, Chromosome, Evaluation>::Individual Individual;
	
protected:
	size_t parentIndex; // used in generateCandidate function
	
public:
	PDe();

protected:
	void performTruncate();
	void performSelection(size_t targetIndex, Individual trial);
	
	Chromosome generateCandidate();
};


template<class Chromosome, class Evaluation>
PDe<Chromosome, Evaluation>::PDe() : parentIndex(0) {
}


template<class Chromosome, class Evaluation>
void PDe<Chromosome, Evaluation>::performSelection(size_t targetIndex, Individual trial) {
	if (trial < this->pop[targetIndex]) 
		this->pop[targetIndex] = trial;
}


template<class Chromosome, class Evaluation>
void PDe<Chromosome, Evaluation>::performTruncate() {
}


template<class Chromosome, class Evaluation>
Chromosome PDe<Chromosome, Evaluation>::generateCandidate() {
	Individual ret;
	if (parentIndex >= this->pop.size())
		parentIndex = 0;
	//this->generateTrialSolution(this->pop, this->indices[parentIndex], ret);
	this->generateTrialSolution(this->pop, parentIndex, ret);
	++parentIndex;
	return ret.chromosome;
}


#endif // DE_H_INCLUDED
