#ifndef DE_H_INCLUDED
#define DE_H_INCLUDED


#include "DeBase.h"
#include "NumericOptimizer.h"

#include <map>


/// ##################################################################################
/// template class De
/// used as an Algorithm policy for NumericOptimizer<...>;
///
template<class Chromosome, class Evaluation>
class De : 
	// inherit base DE definitions
	public DeBase<IndividualStruc<Chromosome, typename Evaluation::Value> >, 
	// make De conform the form of NumericOptimizer
	public NumericOptimizer<De<Chromosome, Evaluation>, Chromosome, Evaluation> 
{
public:
	friend class NumericOptimizer<De<Chromosome, Evaluation>, Chromosome, Evaluation>;
	typedef typename NumericOptimizer<De<Chromosome, Evaluation>, Chromosome, Evaluation>::Individual Individual;
	
public:
	De();

protected:
	void nextGeneration();
	void generationalSelection();
	void performSelection(size_t targetIndex, size_t trialIndex);
};


template<class Chromosome, class Evaluation>
De<Chromosome, Evaluation>::De() {
}


template<class Chromosome, class Evaluation>
void De<Chromosome, Evaluation>::nextGeneration() {
	randomizeIndices(this->indices, this->pop.size());
	
	for (size_t i = 0; i < this->indices.size(); ++i) {
		this->generateTrialSolution(this->pop, this->indices[i], this->children[i]);
		this->eval.normalize(this->children[i].chromosome);
	}
}


template<class Chromosome, class Evaluation>
void De<Chromosome, Evaluation>::generationalSelection() {
	for (size_t i = 0; i < this->indices.size(); ++i) {
		performSelection(this->indices[i], i);
	}
}


template<class Chromosome, class Evaluation>
void De<Chromosome, Evaluation>::performSelection(size_t targetIndex, size_t trialIndex) {
	if (this->children[trialIndex] < this->pop[targetIndex]) 
		this->pop[targetIndex] = this->children[trialIndex];
}


#endif // DE_H_INCLUDED
