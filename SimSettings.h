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

#ifndef SIMSETTINGS_H_INCLUDED
#define SIMSETTINGS_H_INCLUDED

#include <vector>
#include <stdexcept>


/*******************************************************************************************//**
    class SimSettings
    
    all the settings of the simulator and optimizer that are to be read from a file
**/
struct SimSettings {
    /// how should the criteria for the optimization be set
	enum CriteriaMode {
		mode_every_lead_is_criterium = 1,
		mode_leads_sum_is_criterium = 9
	};
	
	/// how is the simulation result to be compared to the measured ECGs
	enum ComparisonMode {
		mode_rms = 1,
		mode_correlation = 2,
		mode_normalization_offset_division_variance = 3,
		mode_vector_correlation = 4
	};
	
	/// base aps that are used to provide base coefficients (k1-k8 for WohlfartPlus)
	std::vector<WohlfartPlus> baseAps;
	/// minimum and maximum values for k (used by the optimizer)
	std::vector<double> kMin, kMax;
	/** how should the unset APs be interpolated
        - "endo-epi"
        - "endo-mid-api"
	**/
	std::string interpolationTypeString;
	/// vector of indices of ks that are to be optimized (other ks are fixed)
	std::vector<int> freeKs;
	
	/// filename of the measured ECGs that are used as targets for the optimization
	std::string optimizationTargetsFname;
	/// comparison mode between simulated and measured ECG 
	/// @see enum ComparisonMode for possible values
	ComparisonMode comparisonMode;
	/// the way criteria are selected
	/// @see enum CriteriaMode for possible values
	CriteriaMode criteriaMode;
	/// additional criteria may be set, similarity between endo and delayed epi
	/// if delay < 0 this criterion is disabled
	double endoEpiMinCriterionDelay;
	/// number of generations in the optimization 
	size_t numGenerations;
	/// size of the population for the optimization
	size_t populationSize;
	/// queue size for the parallel optimization
	int queueSize;
	/// position of the mid layer (in range [0, 1]) between endo and epi
	double midPosition;
	/// additional criterion: for every electrode that is criterion, peak 
	/// positions of simulated and measured ECG are compared (peak = max value)
	bool peakPositionIsCriterion;
	/// calculate fast approximation of an ECG (endo - delayed epi) and use it as an 
	/// extra criterion (compared aginst the first measured ECG)
	bool fastApproxIsCriterion;
	/// the delay of epi used in fast approximation of an ECG
	double fastApproxEpiDelay;
	/// if fast approximation returns result above the limit, regular simulation is not 
	/// executed, only the result of fast approximation is returned instead
	double fastApproxLimit;
	
private:
	template<class El>
	static void outputVector(std::ostream& out, const std::vector<El>& v) {
	    out << "[";
	    if (v.size() > 0) {
            out << v[0];
	    }
	    for (size_t i = 1; i < v.size(); ++i) {
            out << "," << v[i];
        }
	    out << "]";
	}
	
public:
	/// load settings from *.ini file
	void load(const char* fname) {	
	    	
		std::cerr << "***** loading additional simulation parameters *************\n";
		Ini::File ini(fname);
		
		if (ini) {
			std::cerr << " from file: " << fname << "\n";
		} else {
			std::cerr << " file " << fname << " ";
			if (ini.notFound())
				std::cerr << " not found, ";
			if (ini.empty())
				std::cerr << " empty, ";
			if (ini.ok())
				std::cerr << " contains errors.";
			std::cerr << "\n";
		}
		{	
			int apSec = ini.getSectionNumber("wohlfart ap");
			
			/// base ap X x is in [1..inf]
			for (size_t i = 0; ; ++i) {
				std::ostringstream name;
				name << "base ap " << (i+1);
				
				std::vector<double> k;
				Ini::ArrayReader<double, std::vector<double> > reader(k, 9);
				if (ini.loadVar(reader, name.str(), apSec) && (k.size() == 9)) {
					baseAps.push_back(WohlfartPlus());
					baseAps.back().setK(k);
					for (size_t i=0; i<9; ++i)
						std::cerr << " k" << i << "=" << baseAps.back().getK()[i];
					std::cerr << "\n";
				} else
					break;
			}
			if (baseAps.size() < 2)
				throw std::runtime_error("need at least 2 base APs");
			
			/// kMin and kMax
			{
				Ini::ArrayReader<double, std::vector<double> > reader(kMin, 9);
				if (!ini.loadVar(reader, "k min", apSec) || (kMin.size() != 9)) 
					throw std::runtime_error("could not read kMin");
				std::cerr << " k min = ";
				outputVector(std::cerr, kMin);
				std::cerr << "\n";
			}
			{
				Ini::ArrayReader<double, std::vector<double> > reader(kMax, 9);
				if (!ini.loadVar(reader, "k max", apSec) || (kMax.size() != 9)) 
					throw std::runtime_error("could not read kMax");
				std::cerr << " k max = ";
				outputVector(std::cerr, kMax);
				std::cerr << "\n";
			}
			
			/// interpolation type as a string
			/// valid strings are:
			/// 	endo-epi
			/// 	endo-mid-epi
			ini.loadVar(interpolationTypeString, "interpolation", apSec);
			std::cerr << " AP interpolation set to " << interpolationTypeString << " ("
				<< "string validity not checked yet)\n";
			
			/// mid position (only used if mid APs are used)
			/// this variable defines the relative position (between 0[inner border]
			///    and 1[outer border] of AP layers)
			midPosition = 0.5;
			ini.loadVar(midPosition, "mid AP position", apSec);
			if ((midPosition > 1.0) || (midPosition < 0.0)) {
				std::cerr << " warning, mid AP position must be in range [0..1] but is set to "
					<< midPosition << " in settings; using default value of 0.5 instead\n";
				midPosition = 0.5;
			}
			
			/// list of free "Wohlfart k"s (ks that are to be 'optimized')
			{
				Ini::ArrayReader<int, std::vector<int> > reader(freeKs);
				ini.loadVar(reader, "free k", apSec);
				std::cerr << " free Wohlfart koefficients ";
				outputVector(std::cerr, freeKs);
				std::cerr << "\n";
			}
			
			std::cerr << "\n";
		}
		{
			std::cerr << "***** loading optimization parameters **********************\n";
			
			int optSec = ini.getSectionNumber("optimization");
			
			ini.loadVar(optimizationTargetsFname, "targets filename", optSec);
			std::cerr << " targets filename = " << optimizationTargetsFname << "\n";
			
			{
				int crit = 1;
				ini.loadVar(crit, "mode", optSec);
				switch (crit) {
				case mode_every_lead_is_criterium:
				case mode_leads_sum_is_criterium:
					criteriaMode = (CriteriaMode)crit;
					break;
				default:
					throw std::runtime_error("invalid optimization mode");
					break;
				}
				std::cerr << " criterization mode = " << criteriaMode << "\n";
			}
			
			endoEpiMinCriterionDelay = -1;
			ini.loadVar(endoEpiMinCriterionDelay, "endo-epi minimization criterion epi delay", optSec);
			if (endoEpiMinCriterionDelay > 0) 
				std::cerr << "enabling additional criterion - endo-epi minimization with epi " <<
				"delay of " << endoEpiMinCriterionDelay << "\n";
			
			peakPositionIsCriterion = false;
			{
				int temp = 0;
				ini.loadVar(temp, "peak position is criterion", optSec);
				if (temp > 0) 
					std::cerr << "enabling additional criteria - peak position for every base\n";
				peakPositionIsCriterion = temp > 0;
			}
			
			{
				int temp = 0;
				ini.loadVar(temp, "fast approximation is criterion", optSec);
				if (temp > 0) 
					std::cerr << "enabling additional criteron - fast approximation of an ECG\n";
				fastApproxIsCriterion = temp > 0;
				
				fastApproxEpiDelay = 0;
				ini.loadVar(fastApproxEpiDelay, "fast approximation epi delay", optSec);
				
				fastApproxLimit = 2;
				ini.loadVar(fastApproxLimit, "fast approximation limit", optSec);
			}
			
			{
				int cmode = 2;
				ini.loadVar(cmode, "comparison mode", optSec);
				switch(cmode) {
				case mode_correlation:
				case mode_vector_correlation:
				case mode_rms:
				case mode_normalization_offset_division_variance:
					comparisonMode = (ComparisonMode)cmode;
					break;
				default:
					throw std::runtime_error("invalid comparison mode");
					break;
				}
				std::cerr << " ECG comparison mode = " << comparisonMode << "\n";
			}
			
			numGenerations = 100;
			ini.loadVar(numGenerations, "number of generations", optSec);
			std::cerr << " number of generations = " << numGenerations << "\n"; 
			
			populationSize = 0;
			ini.loadVar(populationSize, "population size", optSec);
			std::cerr << " population size = ";
			if (populationSize == 0) std::cerr << "automatic\n";  
			else std::cerr << populationSize << "\n"; 
			
			queueSize = 1;
			ini.loadVar(queueSize, "queue size", optSec);
			std::cerr << " queue size = " << queueSize << "\n";
			
			std::cerr << "\n";
		}
	}
};

#endif // SIMSETTINGS_H_INCLUDED
