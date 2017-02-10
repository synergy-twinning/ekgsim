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
