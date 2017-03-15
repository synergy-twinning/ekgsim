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

#ifndef COLUMNFILEL_H_INCLUDED
#define COLUMNFILEL_H_INCLUDED

/**
    This file contains classes functions for working with .column files
**/

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>


typedef std::vector<double> Column;
typedef std::vector<Column> Columns;


/**
    Column file used for saving vectors of equal length and is readable by matlab. Vectors are saved as columns.
**/
class ColumnFile {
	Columns aps;
	
public:
	ColumnFile() {}
	
	ColumnFile(const char* fname) {
		load(fname);
	}
	
	Columns& columns() {
		return aps;
	}
	
	/// return all columns as a vector of vectors of double
	const Columns& columns() const {
		return aps;
	}
	
	/// returns one column (indexed by p)
	Column& operator[] (size_t p) {
		return aps[p];
	}
	
	const Column& operator[] (size_t p) const {
		return aps[p];
	}
	
	/// number of columns
	size_t numColumns() const {
		return aps.size();
	}
	
	/// returns length of columns (all columns are of equal length)
	size_t columnsLength() const {
		return (aps.size() > 0 ? aps[0].size() : 0);
	}
	
	/// load file
	void load(const char* fname) {
		using namespace std;
		string ext = filenameExtension(fname);
		
		if (ext == ".column") {
			ifstream inFile(fname);
			if (!inFile.is_open()) 
				throw std::runtime_error((string("could not open ") + fname));
				
			// skip comment lines
			while ((inFile.peek() == '#') || (inFile.peek() == '%'))
				inFile.ignore(10000, '\n');
			
			string line;
			getline(inFile, line);
			istringstream lineStr(line);
			vector<double> temp;
			{
				double d;
				while (lineStr >> d) 
					temp.push_back(d);
				lineStr.clear();
			}
			
			if (temp.size() == 0) 
				throw std::runtime_error("File has 0 columns");
			
			aps.resize(temp.size());
			for (size_t i = 0; i < aps.size(); ++i)
				aps[i].resize(0);
			
			while (inFile && lineStr) {
				// if all columns were read successfully, move values to APs
				for (size_t i = 0; i < aps.size(); ++i) {
					aps[i].push_back(temp[i]);
				}
				
				// read line at a time
				getline(inFile, line);
				lineStr.clear();
				lineStr.str(line);
				
				for (size_t i = 0; i < temp.size(); ++i) {
					lineStr >> temp[i];
				}
			}
		} else {
			throw std::runtime_error("ColumnFile::load only works on .column files");
		} 
	}
	
private:
    /// extract extension from filename
	std::string filenameExtension(const std::string& filename) {
		for (int i = (int)filename.size()-1; i >= 0; --i) {
			if (filename[i] == '.') {
				return std::string(filename.begin() + i, filename.end());
			}
		}
		return "";
	}
};


#endif // COLUMNFILEL_H_INCLUDED
