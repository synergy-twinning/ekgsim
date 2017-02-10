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
