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

#include "sim.h"
#include "simlib/sim_lib.h"
#include "nonlinearFit.h"
#include <set>
#include <cmath>
#include "SimSettings.h"


/**
    approximated derivative of AP in the defined time (x)
**/
double der(double x, const SimLib::ActionPotential& ap) {
	//std::cout << "darivative in " << x << " = " << ap(x + 0.1) << " - " << ap(x) << "\n";
	return 100 * (ap(x + 0.01) - ap(x));
}


/*******************************************************************************************//**
    class WohlfartInterpolationEvaluator

    The application of WohlfartInterpolationEvaluator finds a Wohlfart AP params that produce an 
    AP shape, similar to two other given APs (called border APs). The level of similarity is 
    given as a number (called ratio) between 0 and 1, that determines how similar the AP shape
    should be to border AP #1. Similarity to border AP #2 is (1 - ratio). Thus the AP found is 
    in essence an interpolation between the two border APs.
*/
struct WohlfartInterpolationEvaluator {
    /**
        helper class for linear interpolation between two 2D points
    **/
	struct LineConnector {
		/// line parameters (\f$y = kx + n\f$)
		double k, n;
		/// start and end points of the line
		double x1, x2;
		
		/// calculate the value of line in point x
		double operator() (double x) const {
			return k*x + n;
		}
		
		/// transform relative position (0...1) into true position (on x axis)
		double getX(double relativePos) {
			return x1 + (x2 - x1) * relativePos;
		}
		
		/// improve the approximate x towards the true intesection point between the line
		/// and the provided AP; return k of the tangent on the AP in x
		double findIntersection(const SimLib::ActionPotential& ap, double& x) {
			static const double epsY = 0.1;
			static const size_t maxIterations = 10;
			
			// tangent to ap in x
			double l = der(x, ap);
			double m = ap(x) - l * x;
			
			for (size_t cnt = 0; cnt < maxIterations; ++cnt) {
				//std::cout << "[" << x << " " << l << "] ";
				
				// intersection of lines y=k*x+n and y=l*x+m
				if (l != k)
					x = (m-n) / (k-l);
				double y = k*x + n;
				// project the intersection on the AP2
				double yOnAp = ap(x);
				
				// difference between the projected and true intersections
				if (fabs(y - yOnAp) < epsY)
					break;
				
				// tangent in modified x
				l = der(x, ap);
				m = ap(x) - l * x;
			}
			
			if (l != 0)
				return l;
			else 
				return -0.0000000001;
		}
	};
	
	/// @var connectors is a vector of lines connecting points on one AP to points on the other AP
	/// these connectors are used to determine exact points through which an optimzied AP should go
	std::vector<LineConnector> connectors;
	/// @var points is a vector of predefined points over which the optimization performs evaluation
	std::vector<double> points;
	/// @var derivatives is a vector of predefined derivatives over which the optimization performs evaluation
	std::vector<double> derivatives;
	
	/// return sum of squared errors of AP in selected points (including squared errors of differentials)
	/// uses points and \ref derivatives as \ref input
	double operator() (const SimLib::ActionPotential& ap) {
		double sum = 0;
		for (size_t i = 0; i < points.size(); i+=2) {
			sum += sqr(ap(points[i]) - points[i+1]);
		}
		
		for (size_t i = 0; i < derivatives.size(); i+=2) {
			sum += sqr(der(derivatives[i], ap) - derivatives[i+1]);
		}
		
		return sum;
	}
	
	/// set the APs between which the new AP should be searched for
	void setBorderAps(const SimLib::ActionPotential& ap1, const SimLib::ActionPotential& ap2) {
		connectors.clear();
		distributeConnectors(ap1, ap2);
	}
	
	/// For the optimization to work, connectors between the border APs must be formed (border APs
	/// must be set prior to the call of this function). Function tries to distribute connectors
	/// on important points over the whole AP (apd90, highest derivative, ...)
	void distributeConnectors(const SimLib::ActionPotential& ap1, const SimLib::ActionPotential& ap2) {
		// apd_90 will serve as tle last point, while 10ms will serve as the first point
		
		// total number of points to be distributed along the length of the ap
		// must be >= 2
		size_t numPoints = 15;
		
		// find apd_90 of a both aps
		double apd1 = ap1.wohl.apd90();
		double apd2 = ap2.wohl.apd90();
		
		// start x (must be at least 1!)
		size_t startX = 10;
		
		// approximate lengths of both aps (long aps (>700ms) will be truncated)
		// length as in absolute curve length (from start to apd90)
		double len1 = 0.0, len2 = 0.0;
		std::vector<double> ap1y(700), ap2y(700);
		
		static const double xScaleFactor = 0.1;
		{
			for (size_t i = 0; i < 700; ++i) {
				ap1y[i] = ap1(i);
				ap2y[i] = ap2(i);
			}
			
			for (size_t i = startX+1; i < std::min((size_t)700, (size_t)ceil(apd1)); ++i)
				len1 += sqrt(sqr(ap1y[i]-ap1y[i-1]) + xScaleFactor);
			for (size_t i = startX; i < std::min((size_t)700, (size_t)ceil(apd2)); ++i)
				len2 += sqrt(sqr(ap2y[i]-ap2y[i-1]) + xScaleFactor);
		}
		
		double currentLen1 = 0.0, currentLen2 = 0.0;
		size_t i1 = startX;
		size_t i2 = startX;
		LineConnector lc;
		lc.k = 100;
//		std::ofstream pFile("C:/Cygwin/home/mdepolli/points.txt");
		for (size_t i = 0; i < numPoints; ++i) {
			// create a connection between two points, one on each AP, both relatively
			// equally distanced from the AP start (1/numPoints * AP length)
			double targetLen = i * len1 / (numPoints-2);
			for (double l1 = currentLen1; (l1 < targetLen) && (i1 < 700); ++i1)
				l1 += sqrt(sqr(ap1y[i1+1]-ap1y[i1]) + xScaleFactor);
			currentLen1 = targetLen;
			targetLen = i * len2 / (numPoints-2);
			for (double l2 = currentLen2; (l2 < targetLen) && (i2 < 700); ++i2)
				l2 += sqrt(sqr(ap2y[i2+1]-ap2y[i2]) + xScaleFactor);
			currentLen2 = targetLen;
			
			// line throught P (or the points defined with i1 and i2)
			lc.x1 = i1;
			lc.x2 = i2;
			if (i1 == i2)
				lc.x2 += 0.001;
			
			lc.k = (ap1y[i1]-ap2y[i2]) / (lc.x1-lc.x2);
			
			lc.n = ap2y[i2] - lc.k * i2;
//			pFile << lc.x1 << " " << ap1y[i1] << " " << lc.x2 << " " << ap2y[i2] << " " << lc.k << " " << lc.n << "\n";
			
			connectors.push_back(lc);
		}
	}
	
	
	/// setup ratio (a number between 0 and 1, indicating the similarity between the sought for AP and the border AP #1)
	void setupRatio(double ratio) {
		points.clear();
//		std::ofstream pFile("C:/Cygwin/home/mdepolli/points.txt");
		for (size_t i = 0; i < connectors.size(); ++i) {
			points.push_back(connectors[i].getX(ratio));
			points.push_back(connectors[i](points.back()));
//			pFile << points[2*i] << " " << points[2*i+1] << "\n";
		}
	}
	
private:
	/// Create a line connector from two given points (x1, x2) on given APs (ap1, ap2)
	/// y1 and y2 are calculated from x1, x2, and APs. Condition: x1 < x2
	LineConnector connectPoints(const SimLib::ActionPotential& ap1, const SimLib::ActionPotential& ap2, double x1, double x2) {
		static const double epsY = 0.01;
		
		// define a line perpendicular to ap1 at time x as y = k*x + n
		double k = -1.0 / der(x1, ap1);
		double n = ap1(x1) - k * x1;
		
		// *test* define a non perpendicular line
		double y1 = ap1(x1);
		double y2 = ap2(x2);
		k = (y2 - y1) / (x2 - x1);
		n = y1 - k * x1;
		
		// define a tangent on the ap2 as y = l*x + m
		double l = der(x2, ap2);
		double m = ap2(x2) - l * x2;
		
//		std::cout << "connect start " << "y=" << k <<"*x+" << n << "  ";
//		std::cout << "and " << "y=" << l <<"*x+" << m << "\n";
		
		// iteratively find approximate intersection of y=k*x+n and ap2
		for (size_t cnt = 0; cnt < 100; ++cnt) {
			// find the intersection of lines y=k*x+n and y=l*x+m
			double x = (m-n) / (k-l);
			double yOnLine = k*x + n;
			// project the intersection on the AP2
			double yOnAp = ap2(x);
			
			if (fabs(yOnAp - yOnLine) < epsY)
				break;
			
			l = der(x, ap2);
			m = ap2(x) - l * x;
		}
		
		LineConnector line;
		line.k = k;
		line.n = n;
		line.x1 = x1;
		line.x2 = (m-n) / (k-l);
		return line;
	}
	
	/// add a point that AP should go through
	void addPoint(double x, double y) {
		points.push_back(x);
		points.push_back(y);
	}
	
	/// add a derivative at given x, that AP should follow
	void addDerivative(double x, double d) {
		derivatives.push_back(x);
		derivatives.push_back(d);
	}
};


/// *******************************************************************************************//**
/// hidden class SimImplementation (defined in .cpp file)
///
/// class hides EkgSim into .cpp file (for faster compiling)
/// includes actual simulation & evaluation
///
class SimImplementation {
	typedef std::vector<double> Value;
	typedef std::vector<double> Input;
	
	/// sim should not be created at the beginning of ctor (needs some settings first), therefore auto_ptr 
	std::auto_ptr<EkgSim> sim;
	/// target APs
	std::vector<std::vector<double> > targets;
	/// DEPRECATED target offsets (used in comparison function) 
	std::vector<double> targetOffsets;
	
	/// backup of layer aps, used for output
	std::vector<SimLib::ActionPotential> layerAps;
	
	/// approximate simulation result
	std::vector<double> approximateECG;
	/// match between the approximate ECG and the first measurement
	double approxEcgCriteria;
	/// simulationDone is marked false if simulation was not run due to bad result from approximation
	bool simulationDone;
	
	/// types of transmural interpolation
	enum InterpolationType {
		between_endo_epi = 1,
		between_endo_mid_epi,
		interpolation_unknown = 0
	} interpolationType;
	/// set of WohlfartPlus params that are allowed to vary during the optimization
	std::set<char> variableWohlfartParams; 
	
	size_t numDisplacementParams, numWohlfartParams;
	
public:
	SimSettings settings;
	OutputSettings outSettings;
	size_t deducedNumOfCriteria;
	
public:
    /// ctor initializes the simulator, loads its settings, heart model, conduction matrix,
    /// measuring points, and excitation sequence (or simulates it)
	SimImplementation() {
		simulationDone = false;
		/// initialization
		
		/// 1. load and apply simulator settings (one .ini file)
		std::cerr << "***** setting up Ekg Simulator *****************************\n";
		loadSettings("simulator.ini");
		
		/// 2. shape, APs, conduction, measuring points
		sim->loadTransferMatrix();
		sim->loadMeasuringPoints();
		numDisplacementParams = settings.measuringPointsDisplacementIsInput ? sim->numMeasurements()*2 : 0;
		sim->loadShape();
		
		/// 3. either calculate excitation sequence or use the precalculated one
		sim->simExcitationSequence();
		
		/// 5. apply and check settings
		std::cerr << "\n***** applying (and checking) settings *********************\n";
		sim->applySettings();
		if (settings.interpolationTypeString == "endo-epi")
			interpolationType = between_endo_epi;
		else if (settings.interpolationTypeString == "endo-mid-epi")
			interpolationType = between_endo_mid_epi;
		else {
			interpolationType = interpolation_unknown;
			throw std::runtime_error(std::string("unknown interpolation type [") + settings.interpolationTypeString
				+ "]");
		}
		setWohlfartFreeParams(settings.freeKs.begin(), settings.freeKs.end());
		std::cerr << " ok\n";
		
		std::cerr << "\n***** printout of the simulator setup **********************\n";
		sim->printSettings();
		sim->printMessages(std::cout);
		
		/// 5. load targets (targets for the measurements)
		std::cerr << "\n***** loading targets and setting up optimization **********\n";
		loadTargets(settings.optimizationTargetsFname.c_str());
		
		size_t baseDeducedNumOfCriteria = 0;
		switch(settings.criteriaMode) {
		case SimSettings::mode_every_lead_is_criterium:
			baseDeducedNumOfCriteria = std::min(sim->numMeasurements(), targets.size());
			break;
		case SimSettings::mode_leads_sum_is_criterium:
			baseDeducedNumOfCriteria = 1;
			break;
		default:
			break;
		}
		deducedNumOfCriteria = baseDeducedNumOfCriteria;
		if (settings.peakPositionIsCriterion)
			deducedNumOfCriteria += baseDeducedNumOfCriteria;
		if (settings.fastApproxIsCriterion)
			++deducedNumOfCriteria;
		if (settings.endoEpiMinCriterionDelay >= 0)
			++deducedNumOfCriteria;
		std::cerr << " deduced number of criteria = " << deducedNumOfCriteria << "\n";
		std::cerr << "\n";
	}
	
	void loadSettings(const char* fname) {
		// settings for the optimizer
		settings.load(fname);
		
		// settings for the simulator (sadly loading happens twice)
		sim = std::auto_ptr<EkgSim>(new EkgSim);
		sim->loadSettings(fname);
	}
	
	/// evaluate an individual (do a simulation) and return violation
	double eval(const Input& solution, Value& result) {
		static size_t num = 0;
		++num;
		
		std::cout << "\reval " << num << "  ";
		
		PrecisionTimer pt;
		double evalViolation = 0;
		
		{
			ScopeTimer t(pt);
			
			// make sure solution contains enough values to represent all these parameters
			if (solution.size() < numWohlfartParams+numDisplacementParams)
				throw std::runtime_error("Solution does not contain enough values");
			
			if (settings.measuringPointsDisplacementIsInput) {
				// take displacements into a separate vector
				std::vector<EkgSim::PositionVec> displacements(sim->numMeasurements());	
				for (size_t i = 0; i < sim->numMeasurements(); ++i) {
					displacements[i][0] = solution[numWohlfartParams+i*2];
					displacements[i][1] = solution[numWohlfartParams+i*2+1];
					displacements[i][2] = 0;
				}
				
				// use the displacements
				sim->moveMeasuringPoints(displacements);
			}
			
			// now simulate (only wohlfart parameters are needed here, but whole solution can be passed on, since the extra values in it will not effet the solving procedure)
			switch (interpolationType) {
			case between_endo_epi:
				evalViolation = simUsingBorderAps(solution);
				break;
			case between_endo_mid_epi:
				evalViolation = simUsingBorderAndMidAps(solution);
				break;
			default:
				break;
			}
		}
		
//		if (pt.totalSeconds() > 30)
//			std::cout << pt.totalSeconds() << "s " << std::setprecision(26) << solution << "\n";
		evalViolation += calculateFitness(result);
		writeOutputs(solution);
		
		return evalViolation;
	}
	
	template<class Iterator>
	void setWohlfartFreeParams(Iterator begin, Iterator end) {
		std::copy(begin, end, std::inserter(variableWohlfartParams, variableWohlfartParams.begin()));
		
		// determine number of wohlfart parameters
		numWohlfartParams = 0;
		switch (interpolationType) {
		case between_endo_epi:
			numWohlfartParams = variableWohlfartParams.size()*2;
			break;
		case between_endo_mid_epi:
			numWohlfartParams = variableWohlfartParams.size()*3;
			break;
		default:
			break;
		}
	}
	
	/// interpolation type and the free parameters to wohlfart must be set before calling this function
	void getGeneParams(size_t& numGenes, size_t& numCriteria, std::vector<double>& gMin, std::vector<double>& gMax) {
		try {
			switch (interpolationType) {
			case between_endo_epi: {
				gMin.resize(variableWohlfartParams.size()*2);
				gMax.resize(variableWohlfartParams.size()*2);
				size_t gIndex = 0;
				std::cout << "-----------------\n-- gene limits --\n";
				std::cout << " num free params = " << variableWohlfartParams.size() << "\n";
				for (size_t i = 0; i < WohlfartPlus::numParams; ++i) {
					if (variableWohlfartParams.count(i) > 0) {
						gMin[gIndex] = gMin[gIndex + variableWohlfartParams.size()] = settings.kMin[i];
						gMax[gIndex] = gMax[gIndex + variableWohlfartParams.size()] = settings.kMax[i];
						std::cout << " k" << i << "=[" << gMin[gIndex] << "..." << gMax[gIndex] <<
							"]\n";
						++gIndex;
					}
				}
				std::cout << "\n-----------------\n";
				numGenes = variableWohlfartParams.size()*2;
				break;
			}
			case between_endo_mid_epi: {
				gMin.resize(variableWohlfartParams.size()*3);
				gMax.resize(variableWohlfartParams.size()*3);
				size_t gIndex = 0;
				std::cout << "-----------------\n-- gene limits --\n";
				std::cout << " num free params = " << variableWohlfartParams.size() << "\n";
				for (size_t i = 0; i < WohlfartPlus::numParams; ++i) {
					if (variableWohlfartParams.count(i) > 0) {
						gMin[gIndex] = gMin[gIndex + variableWohlfartParams.size()]
							= gMin[gIndex + 2*variableWohlfartParams.size()] = settings.kMin[i];
						gMax[gIndex] = gMax[gIndex + variableWohlfartParams.size()]
							= gMax[gIndex + 2*variableWohlfartParams.size()] = settings.kMax[i];
						std::cout << " k" << i << "=[" << gMin[gIndex] << "..." << gMax[gIndex] <<
							"]\n";
						++gIndex;
					}
				}
				std::cout << "\n-----------------\n";
				numGenes = variableWohlfartParams.size()*3;
				break;
			}
			default:
				break;
			}
			
			numCriteria = deducedNumOfCriteria;
			
			if (settings.measuringPointsDisplacementIsInput) {
				// add 2 coordinates per measuring point as genes 
				numGenes += numDisplacementParams;
				
				if (settings.displacementMin.empty())
					settings.displacementMin.push_back(0);
				if (settings.displacementMax.empty())
					settings.displacementMax.push_back(0);
				// also add value bounds for these genes
				for (int i = 0; i < (int)numDisplacementParams; ++i) {
					gMin.push_back(settings.displacementMin[std::min(i, (int)settings.displacementMin.size()-1)]);
					gMax.push_back(settings.displacementMax[std::min(i, (int)settings.displacementMax.size()-1)]);
				}
			}
		} catch(...) {
			throw std::runtime_error("failed to set gene min and max from the interpolation type and the free parameters list");
		}
	}
	
protected:
    /// get violation of a single gene (how far is its value outside the allowed range)
	double getViolation(SimLib::ActionPotential& ap, size_t kIndex) const {
		if (ap.wohl.getK()[kIndex] < settings.kMin[kIndex])
			return settings.kMin[kIndex] - ap.wohl.getK()[kIndex];
		if (ap.wohl.getK()[kIndex] > settings.kMax[kIndex])
			return ap.wohl.getK()[kIndex] - settings.kMax[kIndex];
		return 0;
	}
	
    /// get violation = sum of violations of genes (deviations from the allowed range of values)
	double getViolation(SimLib::ActionPotential& ap) const {
		double v = 0.0;
		for (size_t i = 0; i < 9; ++i)
			v+= getViolation(ap, i);
		return v;
	}
	
	/// calculates the fitness of the given results
	/// also returns violation (0.0 = no violation, bigger number == stronger violation)
	double calculateFitness(Value& result) {
		double violationSum = 0.0;  // no violation calculation yet here (violation of APs)
									// it is calculated prior to simulation
		
		// evaluate simulation result
		if (settings.criteriaMode == SimSettings::mode_every_lead_is_criterium) {
			if (result.size() == 0)
				result.resize(deducedNumOfCriteria);
			if (result.size() != deducedNumOfCriteria) {
				throw std::runtime_error("error in result.size() - doesn't match deducedNumOfCriteria");
			}
		} else {
			throw std::runtime_error("selected criteria mode is either not implemented yet or invalid");
		}
		
		result[0] = 0;
		
		size_t numMeasurementCriteria = std::min(sim->numMeasurements(), targets.size());
		
		// positions of criteria in result vector:
		// ECG1, [ECG1 peak], ECG2, [ECG2 peak], ..., ECGn, [ECGn peak], 
		//		[fast approx], [endo-epi delay]
		if (settings.fastApproxIsCriterion) {
			result[numMeasurementCriteria] = approxEcgCriteria;
		}
		
		// debug: num simulated measurements and num shape targets
		// std::cout << "numMeas=" << sim->numMeasurements() << ", targets=" << targets.size() << "\n";
		
		// work on single result (sim measurement) at a time
		for (size_t i = 0; i < numMeasurementCriteria; ++i) {
			std::vector<double> simResult;
			ConvolutionResult<double> cr;
			
			if (simulationDone) {			
				// copy result (sim measurements) into local variable
				simResult = sim->getMeasurement(i);
				
				// compare sim measurement to predefined shape
				switch (settings.comparisonMode) {
					case SimSettings::mode_correlation:
						cr = statisticalCorrelationCoeff(simResult, targets[i], 0);
						cr.bestValue = 1.0 - cr.bestValue;
						break;
					case SimSettings::mode_vector_correlation:
						cr = vectorCorrelationCoeff(simResult, targets[i], 0);
						cr.bestValue = 1.0 - cr.bestValue;
						break;
					case SimSettings::mode_rms:
						cr = rms(simResult, targets[i], 0);
						break;
					case SimSettings::mode_normalization_offset_division_variance:
					default:
						cr = devFromLinear(simResult, targets[i], targetOffsets[i], 0);
						break;
				}
			} else {
				simResult = approximateECG;
				cr.bestOffset = 0;
				cr.bestValue = approxEcgCriteria;
			}
			
			if (settings.criteriaMode == SimSettings::mode_every_lead_is_criterium) {
				size_t critIndex = i;
				if (settings.peakPositionIsCriterion) {
					critIndex = 2*i;
					int resultPeak = int(std::max_element(simResult.begin(), simResult.end()) 
						- simResult.begin());
					int targetPeak = int(std::max_element(targets[i].begin(), 
						targets[i].end()) - targets[i].begin());
					result[critIndex+1] = abs(resultPeak - targetPeak);
				}
				result[critIndex] = cr.bestValue;
				// isfinite and isnormal might be GCC only!
				if ((result[critIndex] < 0) || !std::isfinite(result[critIndex]) || !std::isnormal(result[critIndex])) {
					std::cerr << " warning, ConvolutionResult returned wierd value: " 
						<< cr.bestValue << ", causing fitness to be " << result[critIndex] 
						<< "; making correction - setting fitness to a large number\n";
					result[critIndex] = 2e10;
				}
			} else {
				// not yet fully implemented (not of much ues anyways)
				result[0] += cr.bestValue;
			}
		}
		
		// if minimization of the difference between endo and epi is enabled:
		if (settings.endoEpiMinCriterionDelay >= 0) {
			// square difference between the first and last layer APs
			double sd = 0.0;
			for (double t = 1.0; t < 700.0; t += 1.0) {
				double d = layerAps.front()(t + settings.endoEpiMinCriterionDelay) - layerAps.back()(t);
				sd += d*d;
			}
			
			result[result.size() - 1] = sd;
		}
		
		// return violation
		return violationSum;
	}
	
	/// \f$ dest.k = src1.k * (1-ratio) + src2.k * ratio \f$
	void combineAps(SimLib::ActionPotential& dest, const SimLib::ActionPotential& src1, const SimLib::ActionPotential& src2, double ratio) const {
		for (size_t i = 0; i < src1.wohl.numParams; ++i) {
			dest.wohl.getK()[i] = src1.wohl.getK()[i] * (1-ratio) + src2.wohl.getK()[i] * ratio;
		}
	}
	
	/// run approximate simulation and if the result is good enough, full simulation
	void runApproxAndSim() {
		sim->runApproximation(settings.fastApproxEpiDelay, approximateECG);
		
		approxEcgCriteria = 2;
		if (settings.fastApproxIsCriterion || (settings.fastApproxLimit < 2)) {
			ConvolutionResult<double> cr;
			switch (settings.comparisonMode) {
				case SimSettings::mode_correlation:
					cr = statisticalCorrelationCoeff(approximateECG, targets[0], 0);
					cr.bestValue = 1.0 - cr.bestValue;
					break;
				case SimSettings::mode_vector_correlation:
					cr = vectorCorrelationCoeff(approximateECG, targets[0], 0);
					cr.bestValue = 1.0 - cr.bestValue;
					break;
				case SimSettings::mode_rms:
					cr = rms(approximateECG, targets[0], 0);
					break;
				case SimSettings::mode_normalization_offset_division_variance:
				default:
					cr = devFromLinear(approximateECG, targets[0], targetOffsets[0], 0);
					break;
			}
			
			approxEcgCriteria = cr.bestValue;
			//SimLib::exportVector(approximateECG, "approximate.column");
		}
		
		if (approxEcgCriteria > settings.fastApproxLimit) {
			simulationDone = false;
		} else {
			sim->run();
			simulationDone = true;
		}
	}
	
	/// function expects Input to contain (2*number of variable parameters for Wohlfart+ function).
	/// it uses first and last specified base AP for the first and last layer
	double simUsingBorderAps(const Input& solution) {
		assert(settings.baseAps.size() >= 2);
		std::vector<SimLib::ActionPotential> baseAps(sim->requiredAps());
		
		// check if all parameters of the solution will be used and message error 
		// if not (either not all are used, or more are needed)
		
		try {
			//std::cout << " " << std::fixed << solution << "\n";
			size_t solutionI = 0;
			double violation = 0.0;
			
			if (solution.size() != (numDisplacementParams + numWohlfartParams)) {
				std::ostringstream temp;
				temp << "chromosome size does not agree with the combination of the number of free "
					<< "Wohlfart parameters and the selected interpolation procedure (" 
					<< solution.size() << " != " << numDisplacementParams << "+" << numWohlfartParams << ")";
				throw std::runtime_error(temp.str());
			}
			
			SimLib::ActionPotential temp;
			// set first AP
			temp.init(settings.baseAps.front().getK(), 0);
			for (size_t i = 0; i < temp.wohl.numParams; ++i) {
				if (variableWohlfartParams.count(i) > 0) {
					temp.wohl.getK()[i] = solution[solutionI];
					++solutionI;
					violation += getViolation(temp, i);
				}
			}
			baseAps.front().init(temp, 0);
			
			// set second AP			
			temp.init(settings.baseAps.back().getK(), 0);
			for (size_t i = 0; i < temp.wohl.numParams; ++i) {
				if (variableWohlfartParams.count(i) > 0) {
					temp.wohl.getK()[i] = solution[solutionI];
					++solutionI;
					violation += getViolation(temp, i);
				}
			}
			baseAps.back().init(temp, 0);
			
			// set middle APs (interpolate between the border APs)
			WohlfartInterpolationEvaluator wohlEval;
			wohlEval.setBorderAps(baseAps.front(), baseAps.back());
			// experimentally set d:
			double kd[] = {  0,   0,     0, 0.001, 0,   0.00005, 0.0005, 0.01,   0.2};
			SimLib::ActionPotential d;
			d.init(kd, 0);
			for (size_t i = 1; i < baseAps.size()-1; ++i) {
				double ratio = i / double(sim->requiredAps() - 1);
				wohlEval.setupRatio(ratio);
				baseAps[i] = baseAps[i-1];
				combineAps(baseAps[i], baseAps.front(), baseAps.back(), ratio);
				steepestDescend(wohlEval, baseAps[i], d, 0.5, 1e-3, 100);
			}
			
			layerAps = baseAps;
			sim->setApsDestructive(baseAps);
			runApproxAndSim();
			
			return violation;
		} catch (const char* e) {
			throw std::runtime_error(e);
		} catch (std::exception& e) {
			throw;
		} catch (...) {
			throw std::runtime_error("unknown exception occured during simulation");
		}
	}
	
	/// function expects Input to contain (3*number of variable parameters for Wohlfart+ function).
	/// it uses first and last specified base AP for the first and last layer
	double simUsingBorderAndMidAps(const Input& solution) {
		assert(settings.baseAps.size() >= 3);
		size_t numberOfLayers = sim->requiredAps();
		std::vector<SimLib::ActionPotential> baseAps(numberOfLayers);
		
		try {
			size_t solutionI = 0;
			double violation = 0.0;
			size_t absoluteMidPos = (size_t)floor(settings.midPosition * (numberOfLayers-1) + 0.5);
			
			if (solution.size() != (numDisplacementParams + numWohlfartParams)) {
				std::ostringstream temp;
				temp << "chromosome size does not agree with the combination of the number of free "
					<< "Wohlfart parameters and the selected interpolation procedure (" 
					<< solution.size() << " != " << numDisplacementParams << "+" << numWohlfartParams << ")";
				throw std::runtime_error(temp.str());
			}
			
			SimLib::ActionPotential temp;
			// set first AP
			temp.init(settings.baseAps.front().getK(), 0);
			for (size_t i = 0; i < temp.wohl.numParams; ++i) {
				if (variableWohlfartParams.count(i) > 0) {
					temp.wohl.getK()[i] = solution[solutionI];
					++solutionI;
					violation += getViolation(temp, i);
				}
			}
			baseAps.front().init(temp, 0);
			
			// set second (mid) AP
			temp.init(settings.baseAps.back().getK(), 0);
			for (size_t i = 0; i < temp.wohl.numParams; ++i) {
				if (variableWohlfartParams.count(i) > 0) {
					temp.wohl.getK()[i] = solution[solutionI];
					++solutionI;
					violation += getViolation(temp, i);
				}
			}
			baseAps[absoluteMidPos].init(temp, 0);
			
			// set third AP
			temp.init(settings.baseAps.back().getK(), 0);
			for (size_t i = 0; i < temp.wohl.numParams; ++i) {
				if (variableWohlfartParams.count(i) > 0) {
					temp.wohl.getK()[i] = solution[solutionI];
					++solutionI;
					violation += getViolation(temp, i);
				}
			}
			baseAps.back().init(temp, 0);
			
			// set other APs (interpolate between the border and mid APs)
			WohlfartInterpolationEvaluator wohlEval;
			wohlEval.setBorderAps(baseAps.front(), baseAps.back());
			double kd[] = {  0,   0,     0, 0.001, 0,   0.00005, 0.0005, 0.01,   0.2};
			SimLib::ActionPotential d;
			d.init(kd, 0);
			
			for (size_t i = 1; i < numberOfLayers-1; ++i) {
				if (i < absoluteMidPos) {
					double ratio = i / double(absoluteMidPos);
					wohlEval.setBorderAps(baseAps.front(), baseAps[absoluteMidPos]);
					wohlEval.setupRatio(ratio);
					baseAps[i] = baseAps[i-1];
					combineAps(baseAps[i], baseAps.front(), baseAps[absoluteMidPos], ratio);
				} else if (i > absoluteMidPos) {
					double ratio = (i - absoluteMidPos) / 
						double(sim->requiredAps() - absoluteMidPos - 1);
					if ((i - absoluteMidPos) == 1)
						wohlEval.setBorderAps(baseAps[absoluteMidPos], baseAps.back());
					wohlEval.setupRatio(ratio);
					baseAps[i] = baseAps[i-1];
					combineAps(baseAps[i], baseAps[absoluteMidPos], baseAps.back(), ratio);
				}
				if (i != absoluteMidPos)
					steepestDescend(wohlEval, baseAps[i], d, 0.5, 1e-3, 100);
			}
			
			layerAps = baseAps;
			sim->setApsDestructive(baseAps);
			runApproxAndSim();
			
			return violation;
		} catch (const char* e) {
			throw std::runtime_error(e);
		} catch (std::exception& e) {
			throw;
		} catch (...) {
			throw std::runtime_error("unknown exception occured during simulation");
		}
	}
	
    /// when multiple APs are to be saved in a single file, matlab style, set a variable of this type
    /// add all APs to be saved and then save (operator () can also be used, which adds a vector
    /// of apps and saves in one call)
	struct OutputAps {
		std::vector<SimLib::saveVecElement> saveElements;
		std::vector<std::vector<double> > vectors;
		static const size_t maxTime = 700;
		size_t lastIndex;
		
		OutputAps(size_t expectedNum) {
			saveElements.resize(expectedNum);
			vectors.resize(expectedNum);
			lastIndex = 0;
		}
		
		template<class T>
		void addAp(const SimLib::ActionPotential& ap, const T& name) {
			if (lastIndex < vectors.size()) {
				vectors[lastIndex].resize(maxTime);
				saveElements[lastIndex].data = &vectors[lastIndex];
				
				std::ostringstream commentStr;
				commentStr << "k = [" << ap[0];
				for (size_t ki = 1; ki < ap.size(); ++ki)
					commentStr << ", " << ap[ki];
				commentStr << "]";
				saveElements[lastIndex].comment = commentStr.str();
				
				std::ostringstream nameStr;
				nameStr << "ap_" << name;
				saveElements[lastIndex].name = nameStr.str(); 
				
				for (size_t time = 0; time < maxTime; ++time)
					vectors[lastIndex][time] = ap(time);
				++lastIndex;
			}
		}
		
		void operator() (const std::string& fname, const std::vector<SimLib::ActionPotential>& aps, const std::vector<size_t>& indices) {
			for (size_t i = 0; i < indices.size(); ++i) {
				if (indices[i] < aps.size()) {
					addAp(aps[indices[i]], indices[i]);
				}
			}
			saveElements.resize(lastIndex);
			if (saveElements.size() > 0)
				SimLib::exportVectors(saveElements, fname, 0, 1);
		}
		
		void operator() (const std::string& fname, const std::vector<SimLib::ActionPotential>& aps) {
			for (size_t i = 0; i < aps.size(); ++i) {
				addAp(aps[i], i);
			}
			if (saveElements.size() > 0)
				SimLib::exportVectors(saveElements, fname, 0, 1);
		}
	};
	
	/// create files and write the data as specified in output settings (outSettings)
	void writeOutputs(const Input& solution) {
		{
			OutputAps dummy(outSettings.outputCellAps.size());
			dummy("cell_aps.column", sim->getAps(), outSettings.outputCellAps);
		}
		if (outSettings.layerAps) {
			OutputAps dummy(layerAps.size());
			dummy("layer_aps.column", layerAps);
		}
		if (outSettings.result) {
			std::ostringstream s;
			s << "input=" << solution << "; ";
			sim->saveMeasurement(s.str());
		}
	}
	
	
	/// load target curves from .column file 
	void loadTargets(const char* fname) {
		ColumnFile targetFile(fname);
		
		targets.resize(std::max(size_t(1), targetFile.numColumns() - 1));
		targetOffsets.resize(targets.size());
		
		// targets file includes time?
		double targetTimeStep = sim->getSettings().inputActionPotentialsTimeStep; // this seems fishy
							// was inputActionPotentialsTimeStep intendet for this use?
		size_t timeOfs = std::min(size_t(1), targetFile.numColumns() - 1);
		if (timeOfs > 0) {
			targetTimeStep = targetFile[0][1] - targetFile[0][0];
		}
		for (size_t i = timeOfs; i < targetFile.numColumns(); ++i) {
			double tMin, tMax;
			minAndMax(tMin, tMax, targetFile[i]);
			
			// make target range equal 1
			mult(targetFile[i], 1.0 / (tMax - tMin));
			
			// resample to desired timestep and copy to targets
			resample(targetFile[i], targets[i-timeOfs], targetTimeStep / sim->getSettings().simulationTimeStep);
			
			// offset must move target in such a way that its min matches 1 and its max 
			// matches 2
			targetOffsets[i-timeOfs] = 1.0 - (tMin / (tMax - tMin));
		}
		
		if (timeOfs > 0)
			targetTimeStep = targetFile[0][1] - targetFile[0][0];
		
		{
			std::vector<SimLib::saveVecElement>		exportVec(targets.size());
			for (size_t i = 0; i < exportVec.size(); ++i) {
				std::ostringstream		name;
				name << "target " << (i+1);
				exportVec[i].name = name.str();
				exportVec[i].data = &targets[i];
			}
			
			std::string outputDir = ""; // should set this as parameter
			exportVectors(exportVec, (outputDir + "target_chk.column").c_str(), targetFile[0][0], 
				sim->getSettings().simulationTimeStep);
		}
		//timeMultiplier = sim.getSettings().simulationTimeStep;
	}
};


void testWohlInterpolation() {
	{
		SimLib::ActionPotential ap1, ap2, ap12, d;
		WohlfartInterpolationEvaluator wohlEval;
		
		double k1[] = {-86.0, 2.5, 100, 0.9,   0.1, 0.002,   0.03,   1,    325};
		double k2[] = {-86.0, 2.5, 100, 0.85,  0.1, 0.003,   0.033,  1.2,  355};
		double kd[] = {  0,   0,     0, 0.001, 0,   0.00005, 0.0005, 0.01,   0.2};
		
		ap1.init(k1, 0);
		ap2.init(k2, 0);
		d.init(kd, 0);
		ap12.init(ap1, 0);
		wohlEval.setBorderAps(ap1, ap2);
		
		std::ofstream kFile("C:/Cygwin/home/mdepolli/k.txt");
		for (size_t i=0; i<9; ++i)
			kFile << ap1.getK()[i] << "\t";
		kFile << "\n";
		
		for (size_t i = 1; i < 10; ++i) {
			wohlEval.setupRatio(0.1 * i);
			steepestDescend(wohlEval, ap12, d, 0.5, 1e-3, 100);
			for (size_t i=0; i<9; ++i)
				kFile << ap12.getK()[i] << "\t";
			kFile << "\n";
		}
		
		for (size_t i=0; i<9; ++i)
			kFile << ap2.getK()[i] << "\t";
		kFile << "\n";
	}
}


// *******************************************************************************************//**
// struct OptimizationFunction
//
OptimizationFunction::OptimizationFunction(size_t propertiesL) : impl(new SimImplementation) {
	valueLen = impl->deducedNumOfCriteria;
	propertiesLen = propertiesL;
	geneBounds = false;
}


OptimizationFunction::~OptimizationFunction() {
	delete impl;
}


void OptimizationFunction::setup(const OutputSettings& outSet) {
	impl->outSettings = outSet;
}


void OptimizationFunction::operator() (const Input& solution, Value& result, double& violation, Properties& properties) const  {
	if (valueLen > 0) TypeWrapper::resize(result, valueLen);
	if (propertiesLen > 0) TypeWrapper::resize(properties, propertiesLen);
	
	violation = impl->eval(solution, result);
}


void OptimizationFunction::getGeneParams(size_t& numGenes, size_t& numCriteria, std::vector<double>& gMin, std::vector<double>& gMax) {
	impl->getGeneParams(numGenes, numCriteria, gMin, gMax);
}


void OptimizationFunction::getEvolutionaryParams(size_t& numGenerations, size_t& popSize, int& queueSize) {
	numGenerations = impl->settings.numGenerations;
	popSize = impl->settings.populationSize;
	queueSize = impl->settings.queueSize;
}
