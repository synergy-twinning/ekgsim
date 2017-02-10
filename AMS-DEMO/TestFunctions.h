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
