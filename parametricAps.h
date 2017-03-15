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

#ifndef PARAMETRICAPS_H_INCLUDED
#define PARAMETRICAPS_H_INCLUDED


/// *******************************************************************************************//**
/// struct ParametricAPs
/// 
/// Structure which holds 3 Wohlfart APs and offers different possible interpolations into full 
/// range of transmural APs (APs across the myocardium wall)
/// this class is DEPRICATED, should be used instead
/// 
struct ParametricAPs {
	WohlfartPlus endo, mid, epi;
	size_t numAps;
	
	ParametricAPs(size_t num = 12) {
		numAps = num;
		defaultFactors();
	}
	
	/// generate vector of APs by interpolating from endo and epi APs
	void generateApsVectEndoEpiInterpolation(std::vector<std::vector<double> >& ret, std::vector<double>& midKs) {
		// calculate values for k that will produce function resembling linear interpolation
		// between endo and epi functions
		interpolateK(midKs);
		
		ret.resize(numAps + 1);
		for (size_t i = 0; i < ret.size(); ++i) {
			ret[i].reserve(700);
		}
		
		std::vector<WohlfartPlus> w;
		w.resize(numAps - 2);
		double kTemp[9];
		std::copy(endo.getK(), endo.getK() + 9, kTemp);
		
		for (size_t i = 0; i < w.size(); ++i) {
			kTemp[4] = midKs[4*i];
			kTemp[5] = midKs[4*i + 1];
			kTemp[6] = midKs[4*i + 2];
			kTemp[7] = midKs[4*i + 3];
			w[i].setK(kTemp);
		}
		
		for (double time = 0; time < 700.0; time += 1.0) {
			ret[0].push_back(time);
			ret[1].push_back(endo[time]);
			ret[numAps].push_back(epi[time]);
			for (size_t i = 2; i < numAps; ++i) {
				// get value of AP[i]
				ret[i].push_back(w[i-2][time]);
			}
		}
	}
	
	/// generate vector using only endo and epi APs
	void generateApsVectMidInterpolation(std::vector<std::vector<double> >& ret) {
		// calculate values for k that will produce function resembling linear interpolation
		// between mid, and endo or epi functions
		
		// create mid mixed APs
		std::vector<double> ks;
		interpolateKwMid(ks);
		
		ret.resize(numAps + 1);
		for (size_t i = 0; i < ret.size(); ++i) {
			ret[i].reserve(700);
		}
		
		std::vector<WohlfartPlus> w;
		w.resize(numAps);
		double kTemp[9];
		std::copy(endo.getK(), endo.getK() + 9, kTemp);
		
		for (size_t i = 0; i < numAps; ++i) {
			kTemp[4] = ks[4*i];
			kTemp[5] = ks[4*i + 1];
			kTemp[6] = ks[4*i + 2];
			kTemp[7] = ks[4*i + 3];
			w[i].setK(kTemp);
		}
		
		for (double time = 0; time < 700.0; time += 1.0) {
			ret[0].push_back(time);
			for (size_t i = 0; i < numAps; ++i) {
				// get value of AP[i]
				ret[i+1].push_back(w[i][time]);
			}
		}
	}
	
	/// generate time vector and 12 AP vectors of length 700ms (0-699) in steps of 1ms
	/// mid APs are calculated through time interpolation
	void generateApsVect(std::vector<std::vector<double> >& ret) {
		std::vector<double> timeFactors;
		timeFactors.resize(numAps);
		double value;
		
		double mid90 = mid.apd90();
		double endo90 = endo.apd90();
		double epi90 = epi.apd90();
		for (size_t i = 0; i < numAps; ++i) {
			if (factors[i] < 0) {
				timeFactors[i] = -1;
			} else {
				if (i < midPosition) {
					timeFactors[i] = mid90 / (endo90 + factors[i] * (mid90 - endo90));
				} else {
					timeFactors[i] = mid90 / (epi90 + factors[i] * (mid90 - epi90));
				}
			}
		}
		
		ret.resize(numAps + 1);
		for (size_t i = 0; i < ret.size(); ++i) {
			ret[i].reserve(700);
		}
		
		for (double time = 0; time < 700.0; time += 1.0) {
			ret[0].push_back(time);
			for (size_t i = 0; i < numAps; ++i) {
				// get value of AP[i]				
				if (timeFactors[i] < 0) {
					if (i < midPosition) 
						value = endo[time];
					else 
						value = epi[time];
				} else {
					value = mid[time * timeFactors[i]];
				}
				
				ret[i+1].push_back(value);
			}
		}
	}
	
	void outputWohlfart(std::ostream& out, WohlfartPlus& w) {
		for (size_t i = 0; i < 7; ++i)
			out << w.getK()[i] << ", ";
		out << w.getK()[7];
	}
	
	void generateFile(const char* fname = "generatedAP.column") {
		std::vector<std::vector<double> > vect;
		generateApsVect(vect);
		
		double mid90 = mid.apd90();
		double endo90 = endo.apd90();
		double epi90 = epi.apd90();
		std::vector<double> timeFactors;
		timeFactors.resize(numAps);
		
		for (size_t i = 0; i < numAps; ++i) {
			if (factors[i] < 0) {
				timeFactors[i] = -1;
			} else {
				if (i < midPosition) {
					timeFactors[i] = mid90 / (endo90 + factors[i] * (mid90 - endo90));
				} else {
					timeFactors[i] = mid90 / (epi90 + factors[i] * (mid90 - epi90));
				}
			}
		}
		
		std::ofstream file(fname);
		file << "# EkgSim2 generated file\n";
		file << "# k_endo = ";
		outputWohlfart(file, endo);
		file << "\n";
		file << "# k_mid  = ";
		outputWohlfart(file, mid);
		file << "\n";
		file << "# k_epi  = ";
		outputWohlfart(file, epi);
		file << "\n";
		file << "# apd90 = " << endo90 << " (endo), " << mid90 << " (mid), " << epi90 << " (epi)\n";
		file << "# mids are time-scaled: ";
		
		for (size_t i = 1; i < numAps-1; ++i) {
			file << timeFactors[i] << "  ";
		}
		file << "\n";
		for (size_t i = 0; i < vect[0].size(); ++i) {
			file << vect[0][i];
			for (size_t a = 1; a < vect.size(); ++a) {
				file << " " << vect[i][a];
			}
			file << "\n";
		}
	}
	
	void generateFileEpiEndoInterpolation(const char* fname = "generatedAP.column") {
		std::vector<std::vector<double> > dd;
		std::vector<double> ks;
		generateApsVectEndoEpiInterpolation(dd, ks);
		
		double endo90 = endo.apd90();
		double epi90 = epi.apd90();
		std::vector<double> mid90(numAps - 2);
		
		std::ofstream file(fname);
		file << "# EkgSim2 generated file\n";
		file << "# k_endo = ";
		outputWohlfart(file, endo);
		file << "\n";
		
		for (size_t i = 0; i < mid90.size(); ++i) {
			WohlfartPlus m = endo;
			m.getK()[4] = ks[i*4 + 0];
			m.getK()[5] = ks[i*4 + 1];
			m.getK()[6] = ks[i*4 + 2];
			m.getK()[7] = ks[i*4 + 3];
			mid90[i] = m.apd90();
			file << "# k_mid_" << (i+1) << " = ";
			outputWohlfart(file, m);
			file << "\n";
		}
		
		file << "# k_epi  = ";
		outputWohlfart(file, epi);
		file << "\n";
		
		file << "# apd90 = " << endo90 << " (endo)";
		for (size_t i = 0; i < mid90.size(); ++i)
			file << ", " << mid90[i];
		file << " (mids), " << epi90 << " (epi)\n";
		file << "# - mids are interpolated between endo and epi";
		file << "\n";
		for (size_t time = 0; time < 600.0; ++time) {
			for (size_t i = 0; i < numAps + 1; ++i) {
				file << " " << dd[i][time];
			}
			file << "\n";
		}
	}
	
	void generateFileMidInterpolation(const char* fname = "generatedAP.column") {
		std::vector<std::vector<double> > dd;
		generateApsVectMidInterpolation(dd);
		
		std::ofstream file(fname);
		file << "# EkgSim2 generated file\n";
		file << "# k_endo = ";
		outputWohlfart(file, endo);
		file << "\n";
		file << "# k_mid  = ";
		outputWohlfart(file, mid);
		file << "\n";
		file << "# k_epi  = ";
		outputWohlfart(file, epi);
		file << "\n";
		double endo90 = endo.apd90();
		double epi90 = epi.apd90();
		double mid90 = mid.apd90();
		file << "# apd90 = " << endo90 << " (endo), " << mid90 << " (mid), " << epi90 << " (epi)\n";
		file << "# mids are interpolated towards endo and epi \n";
		for (size_t time = 0; time < 600.0; ++time) {
			for (size_t i = 0; i < numAps + 1; ++i) {
				file << " " << dd[i][time];
			}
			file << "\n";
		}
	}
	
	void setFactors(std::vector<double>& f, int pos) {
		factors.resize(f.size());
		midPosition = pos;
		std::copy(f.begin(), f.end(), factors.begin());
	}
	
private:
	// scaling factors
	std::vector<double> factors;
	// position of mid AP
	size_t midPosition;
	
	void defaultFactors() {
		midPosition = 4;
		factors.resize(numAps);
		factors[0] = factors[numAps-1] = -1;
		for (size_t i = 1; i < numAps-1; ++i)
			factors[i] = 1.0;
		/*
		double eck[12] = {
							-1.0,    0.78571, 0.85714, 0.92857,
							1.00000, 0.88889, 0.77778, 0.61111,
							0.50000, 0.38889, 0.27778, -1.0};
		std::copy(eck, eck+numAps, factors.begin());
		*/
	}
	
	/// using endo and epi, produce mid APs (k5-k8 only)
	/// DEPRICATED
	void interpolateK(std::vector<double>& midK) {
		midK.resize((numAps-2)*4);
		
		// calculate minimum diferential for each of the two input APs
		double minDifEndo = (endo[endo.getK()[7]-1] - endo[endo.getK()[7]+1]) * 0.5;
		double minDifEpi = (epi[epi.getK()[7]-1] - epi[epi.getK()[7]+1]) * 0.5;
		
		for (size_t i = 0; i < numAps-2; ++i) {
			WohlfartPlus temp = interpolateBetween(endo, minDifEndo, epi, minDifEpi, (double(numAps-2-i) / (numAps-1)));
			
			midK[i*4 + 0] = temp.getK()[4];
			midK[i*4 + 1] = temp.getK()[5];
			midK[i*4 + 2] = temp.getK()[6];
			midK[i*4 + 3] = temp.getK()[7];
		}
	}
	
	/// using endo and epi and mid, produce other transmural APs (uses k5-k8 only)
	/// DEPRICATED
	// midK[0,1,2,...] = k5, k6, k7 (all first mid), ...
	void interpolateKwMid(std::vector<double>& midK) {
		midK.resize(numAps*4);
		
		// calculate minimum diferential of input APs
		double minDifEndo = (endo[endo.getK()[7]-1] - endo[endo.getK()[7]+1]) * 0.5;
		double minDifEpi = (epi[epi.getK()[7]-1] - epi[epi.getK()[7]+1]) * 0.5;
		double minDifMid = (mid[mid.getK()[7]-1] - mid[mid.getK()[7]+1]) * 0.5;
		
		for (size_t i = 0; i < numAps; ++i) {
			WohlfartPlus temp = (factors[i] <= 0 ?
				(i < midPosition ? endo : epi) :
				(i < midPosition ? 
					interpolateBetween(mid, minDifMid, endo, minDifEndo, factors[i])
				: 
					interpolateBetween(mid, minDifMid, epi, minDifEpi, factors[i])
				)
			);
			
			midK[i*4 + 0] = temp.getK()[4];
			midK[i*4 + 1] = temp.getK()[5];
			midK[i*4 + 2] = temp.getK()[6];
			midK[i*4 + 3] = temp.getK()[7];
		}
	}
	
	/// takes two WohlfartPlus APs (a and b) and make fast jand inaccurate intepolation between them on the distance ka = (0...1) from AP a and 1-ka from AP b
	/// minimum differentials must be profided for both APs
	/// DEPRICATED
	WohlfartPlus interpolateBetween(const WohlfartPlus& a, double minDifA, const WohlfartPlus& b, double minDifB, double ka) {
		WohlfartPlus ab = a;
		
		double kb = 1.0 - ka;
		ab.getK()[4] = (a.getK()[4]*ka + b.getK()[4]*kb);
		if (minDifA != minDifB) {
			// linear interpolation of k6 creates linear interpolation of min linear coefficient
			// of AP derivative when linear interpolation of angle (atan(k)) is actually needed
			double targetMinDif = tan(atan(minDifA)*ka + atan(minDifB)*kb);
			ab.getK()[6] = (a.getK()[6] + (targetMinDif - minDifA)/(minDifB - minDifA) * (b.getK()[6] - a.getK()[6]));
			ab.getK()[5] = ab.getK()[6] / (ka*a.getK()[6]/a.getK()[5] + kb*b.getK()[6]/b.getK()[5]);
		} else {
			ab.getK()[6] = a.getK()[6];
			ab.getK()[5] = a.getK()[5]*ka + b.getK()[5]*kb;
		}
		ab.getK()[7] = a.getK()[7]*ka + b.getK()[7]*kb;
		
		return ab;
	}
};




#endif // PARAMETRICAPS_H_INCLUDED
