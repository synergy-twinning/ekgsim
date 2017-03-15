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

#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED


#include "support.h"
#include <stdexcept>


namespace SimLib {
	
	// ***********************************************************************************************
	// struct ShapeElement
	//
	/// Holds data for one element (or cell) of shape model. Data consists of cell leyer, AP, and 
	/// excitation delay (time lag between the start of simulation and the given cell excitation time)
	/// Struct also contains functions for loading element data from file
	///
	struct ShapeElement {
		size_t								layer;
		double								excitationDelay;
		size_t								apIndex;
		
		static const size_t					layerStartingPoint = 0x1000;
		
		ShapeElement(size_t m = 0) : layer(m), excitationDelay(0), apIndex(m) {}
		
		friend void operator>> (std::istream& in, ShapeElement& el) {
			in >> el.layer;
		}
		
		friend void loadShape(std::istream& in, ShapeElement& el) {
			in >> el.layer;
		}
		
		friend void loadDelay(std::istream& in, ShapeElement& el) {
			in >> el.excitationDelay;
		}
	};


    /// Shape matrix is the data structure that holds shape model data as a 3D matrix of ShapeElements
	typedef Pattern::Hypermatrix<ShapeElement, 3, std::vector<ShapeElement> > ShapeMatrix;
	/// Excitation Transfer Matrix holds data of excitation velocities (speeds of excitation transfer between different layers and inside a layer)
	typedef Pattern::Hypermatrix<double, 2, std::vector<double> > ExcitationTransferMatrix;
	/// MatrixIndex is a data type used for indexing 3D matrices
	typedef Pattern::Vector<size_t, 3> MatrixIndex;

    /// save data from ShapeMatrix in a file (*.matrix), takes additional parameter that tells it weather layers (true) or excitation sequence should be saved (false)
	void exportShape(const ShapeMatrix& shape, const string& filename, bool excit = false);


    /// helper funciton for shapeLoad
	template<class T>
	void defaultLoad(std::istream& in, T& el) {
		in >> el;
	}


    /// function for loading ShapeMatrix (loads layers and excitation start)
    template<class Ch, class ChTr>
	void shapeLoad(std::basic_istream<Ch, ChTr>& in, ShapeElement& el) {
		int layer;
		in >> layer;
		el.layer = (layer < 0 ? ShapeElement::layerStartingPoint - layer : (size_t)layer);
	}


    /// helper function for shapeLoad
	template<class Mtx>
	void importMatrix(Mtx& matrix, const string& filename, void loadFunc(std::istream&, typename Mtx::DataType&) = &SimLib::defaultLoad) {
		string								ext = filenameExtension(filename);
				
		if (ext == ".matrix") {
			ifstream							inFile(filename.c_str());
			
			if (!inFile.is_open()) {
				throw std::runtime_error(string("could not open file ") + filename);
			}
			
			typedef typename Mtx::SizeVector SizeVector;
			typedef typename Mtx::DataType DataType;
			
			// skip comment line
			inFile.ignore(100000000, '\n');
			SizeVector						newSize;
			for (size_t i = 0; i < SizeVector::dimension; ++i)
				newSize[i] = 1;
				
			// read matrix header (dimensionality)
			string							dimStr;
			size_t								dimension;
			inFile >> dimStr;
			if ( ((dimStr[1] == 'D') || (dimStr[1] == 'd')) && ((dimStr[0] > '0') && 
				(dimStr[0] <= ('0' + SizeVector::dimension))) ) {
				// dimension ok
				dimension = dimStr[0] - '0';
			} else {
				throw std::runtime_error("error while loading matrix, could not read dimensionality");
			}
			
			// read matrix header (size vector)
			int writeDim = matrix.dimensionality() - 1;
			newSize[0] = newSize[1] = 1;
			inFile >> newSize[writeDim];
			for (--dimension; dimension > 0; --dimension) {
				inFile.ignore(100, 'x');
				--writeDim;
				inFile >> newSize[writeDim];
			}
			// read matrix header (separator type)
			string							separatorTypeStr;
			getline(inFile, separatorTypeStr, '\n');
			istringstream					separatorTypeSS(separatorTypeStr);
			string							separatorType;
			separatorTypeSS >> separatorType;
			char								separator = (separatorType == "tab" ? '\t' : (
				(separatorType == "comma") || (separatorType == ",") ? ',' : 
				(separatorType == "space" ? ' ' : 0)));
			
			if (inFile.fail()) {				
				throw std::runtime_error(string("reading file ") + filename + " failed (while reading size)");
			}
			matrix.resize(newSize);
			
			// read data
			SizeVector						counter;
			for (size_t i = 0; i < SizeVector::dimension; ++i)
				counter[i] = 0;
			bool								okToLoop = true;
			do {
				char							ch;
				DataType						d;
				//std::cout << counter[2] << "," << counter[1] << "," << counter[0] << "\n"; // DEBUG: see coordinates for loading
				
				switch(separator) {
				case ',':
				case '\t':
				case ' ':
					loadFunc(inFile, d);
					matrix[counter] = d;
					break;
				default:				
					inFile >> ch;
					
					if ((ch == 'x') || (ch == 'X')) {
						matrix[counter] = ShapeElement::layerStartingPoint + 1;
					} else if (ch >= 'a') {
						matrix[counter] = ch - 'a' + 10;
					} else if (ch >= 'A') {
						matrix[counter] = ch - 'A' + 10;
					} else {
						matrix[counter] = ch - '0';
					}
					break;
				}
				
				size_t							incLevel = 0;
				for (size_t targetDim = SizeVector::dimension; targetDim > 0; --targetDim) {
					++counter[targetDim-1];
					if (counter[targetDim-1] >= newSize[targetDim-1]) {
						++incLevel;
						counter[targetDim-1] = 0;
						if (targetDim == 1)
							okToLoop = false;
					} else
						break;
				}
				
				switch(separator) {
				case ',':
					if (incLevel == 0)
						inFile.ignore(100, ',');
				default:	
					break;
				}
			} while (okToLoop);
			
			if (inFile.fail()) {
				throw std::runtime_error(string("reading file ") + filename + " failed (while reading data)");
			}
		} else {
			throw std::runtime_error("can only import from .matrix files for the time being");
		}
	}
	
}


#endif // MATRIX_H_INCLUDED
