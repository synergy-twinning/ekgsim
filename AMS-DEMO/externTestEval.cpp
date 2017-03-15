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

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <iostream>


template<class T>
T sqr(T a) {
	return a*a;
}


template<class Function> 
class FunctionFit {
	std::vector<double> targetValues;
	std::vector<double> parameters;
	std::vector<double> result;
	size_t pairSize;
	std::string directory;
	
public:
	void loadTarget(const char* fname) {
		std::ifstream binFile((fname + std::string(".bin")).c_str(), std::ios::binary);
		if (binFile) {
			// load from binary file
			size_t tempSize;
			binFile.read((char*)&pairSize, 4);
			binFile.read((char*)&tempSize, 4);
			targetValues.resize(tempSize);
			
			for (size_t i = 0; i < targetValues.size(); ++i)
				binFile.read((char*)&targetValues[i], sizeof(double));
			
		} else {
			// load from text file
			pairSize = 2;
			std::ifstream txtFile(fname);
			double x, y;
			while (txtFile) {
				txtFile >> x;
				txtFile.ignore(1);
				txtFile >> y;
				if (!txtFile)
					break;
				txtFile.ignore(1);
				targetValues.push_back(x);
				targetValues.push_back(y);
			}
			
			std::ofstream binOutFile((fname + std::string(".bin")).c_str(), std::ios::binary);
			binOutFile.write((char*)&pairSize, 4);
			size_t tempSize = targetValues.size();
			binOutFile.write((char*)&tempSize, 4);
			for (size_t i = 0; i < targetValues.size(); ++i)
				binOutFile.write((char*)&targetValues[i], sizeof(double));
		}
	}
	
	void loadInput() {
		std::string fname("input.txt");
		if (directory != "")
			if ((directory[directory.size()-1] == '\\') || (directory[directory.size()-1] == '/'))
				fname = directory + fname;
			else 
				fname = directory + "/" + fname;
		std::ifstream file(fname.c_str());
		parameters.reserve(100);
		for(;;) {
			if (file.peek() == '#') 
				file.ignore(1000, '\n');
			else {
				double temp;
				file >> temp;
				if (!file)
					break;
				else
					parameters.push_back(temp);
			}
		}
	}
	
	void writeOutput() {
		std::string fname("output.txt");
		if (directory != "")
			if ((directory[directory.size()-1] == '\\') || (directory[directory.size()-1] == '/'))
				fname = directory + fname;
			else 
				fname = directory + "/" + fname;
		std::ofstream file(fname.c_str());
		
		file << "# writen by externTestEval\n";
		for (size_t i = 0; i < result.size(); ++i) {
			file << result[i] << " ";
		}
	}
	
	void eval(const std::string& dir) {
		directory = dir;
		loadInput();
		Function func;
		result.clear();
		result.resize(2, 0.0);
		
		for (size_t i = 0; i < targetValues.size(); i += 2) {
			result[0] += (sqr(func(parameters, targetValues[i]) - targetValues[i+1]) / targetValues[i+1]);
			result[1] += func(parameters, targetValues[i]) - targetValues[i+1];
		}
		
		writeOutput();
	}
};


struct TestFunc {
	double operator() (const std::vector<double>& params, double in) {
		return params[0] + params[1]*pow(M_E, params[2] * in);
	}
};


int main(int argc, char ** argv) {
	FunctionFit<TestFunc> fit;
	std::string dir("");
	if (argc > 1)
		dir = argv[1];
	fit.loadTarget("target.txt");
	fit.eval(dir);
}
