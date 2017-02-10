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
