#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED


#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cassert>
#include <stdexcept>
#include "Array.h"


class IndividualsFile {
	std::string fileName;
	
	enum {
		read_mode_chromosome = 1,
		read_mode_violation = 2,
		read_mode_criteria = 4,
		read_mode_properties = 8
	};
	
public:
	struct AnalysisResult {
		struct Generation {
			int populationSize;
			int cardinalNum;
			bool errors;
			
			Generation(int c = -1, int s = 0) : 
				populationSize(s), cardinalNum(c), errors(false) 
			{
			}
		};
		
		std::vector<Generation> generations;
		bool orderedGenerations;
		bool staticPopulationSize;
		bool errors;
		bool bigError;
		bool failedToOpenFile;
		int fileMode;
		
		AnalysisResult() : 
			orderedGenerations(true), staticPopulationSize(true), errors(false),
			bigError(false), failedToOpenFile(false), fileMode(-1)
		{
		}
	};
	
	template<class Ind>
	struct FileLine {
		Ind individual;
		int rank;
		double evalTime, lifetime;
	};
	
public:
	IndividualsFile(const char* fname) : fileName(fname) {}
	
	// although Ind is a template parameter, some restrictions apply
	// it is assumed that storage classes for chromosome, criteria and properties
	// are vector types (have members clear(), push_back(x), operator[i])
	template<class Ind, class Stream>
	int readLine(FileLine<Ind>& fline, Stream& stream, int minGen, int& readMode) {
		// first number is always generation, no matter read mode
		int g;
		stream >> g;
		if (!stream) 
			return ~0;
		
		// filter unwanted generations
		if (g < minGen) {
			stream.ignore(100000, '\n');
			return g;
		}
		
		// read chromosome
		if ((readMode & read_mode_chromosome) > 0) {
			size_t read = readField(fline.individual.chromosome, stream);
			if ((read == 0) || !stream) {
				std::cerr << "IndividualsFile::readLine error: " << read << " elements of chromosome read, stream status = " << (bool)stream << "\n";
				return -g;
			}
		}
		
		// violation
		if ((readMode & read_mode_violation) > 0) {
			stream >> fline.individual.violation;
			if (!stream)
				std::cerr << "IndividualsFile::readLine error: " << "unable to read violation\n";
		}
		
		// properties
		if ((readMode & read_mode_properties) > 0) {
			size_t read = readField(fline.individual.getProperties(), stream);
			if (read < fline.individual.getProperties().size())
				fline.individual.getProperties().resize(read);
			if (!stream) {
				std::cerr << "IndividualsFile::readLine error: " << read << " elements of properties read, stream status = " << (bool)stream << "\n";
				return -g;
			}
		}
		
		// criteria
		if ((readMode & read_mode_criteria) > 0) {
			size_t read = readField(fline.individual.criteria, stream);
			if ((read == 0) || !stream) {
				// if could not read criteria then perhaps file is of older type, not 
				// including properties; swap properties and criteria and continue
				
				stream.clear();
				std::cerr << "IndividualsFile::readLine warning: could not read individual in full"
					<< ", one field missing (assuming properties are the missing field)\n";
                
				readMode &= ~read_mode_properties;
				//return -g;
			}
		}
		
		// evaluations file
		if (stream) {
			std::string theRestString;
			std::getline(stream, theRestString);
			std::istringstream theRest(theRestString);
			
			fline.rank = -1;
			fline.evalTime = fline.lifetime = 0.0;
			theRest >> fline.rank >> fline.evalTime >> fline.lifetime;
			if (!stream) {
				stream.clear();
			}
		}
		
		// skip to the end of the line (NEW: done in previous block)
		//stream.ignore(1000, '\n');
		
		if (stream)
			return g;
		else
			return -g;
	}
	
	
	// read single field from stream (terminates at line end or field separator)
	// Field must be indexable type with pre-reserved memory to hold enough elements;
	// stream is read until a delimiter character is found (field size is not checked);
	// returns number of numbers read
	template<class Field, class Stream>
	size_t readField(Field& field, Stream& stream) {
		field.clear();
		for (size_t i = 0; ;) {
			// field may begin with a number or whitespace followed by a number, so try
			// reading that first
			double num;
			stream >> num;
			if (!stream) {
				// failed to read a number				
				stream.clear();
				
				if (i == 0) {
					// this is the beginning of the field, see if it is a delimiter
					// (if not -> error)
					char delimiter;
					stream.get(delimiter);
					// whitespace and vector start delimiters are not an error
					if ((delimiter == ' ') || (delimiter == '\t') || 
						(delimiter == '<') || (delimiter == '(') || (delimiter == '[')) 
						continue;
					// delimiter can also be a slash, marking an empty entry
					if ((delimiter == '/'))
						continue;
					// another option for empty entry
					if ((delimiter == '>') || (delimiter == ')') || (delimiter == ']'))
						return i;
					
					// nothing of the above? signal nothing has been read
					return 0;
				} else {
					// not the first number of the field -> unrecoverable error
					return i;
				}
				continue;
			} else {
				// successful read
				field.push_back(num);
				++i;
			}
			// read dividing character (space, comma, semicolon, ...) or delimiter 
			// (tab, newline, right prantheses ...)
			char delimiter;
			stream.get(delimiter);
			if ((delimiter == '\t') || (delimiter == '\n') || 
				(delimiter == '>') || (delimiter == ')') || (delimiter == ']')) {
				// end of field
				if (delimiter == '\n')
					stream.unget();
				return i;
			}
		}
		// this part of code should be unreachable ...
		return 0;
	}
	
	// returns false if file not found or not right format, true otherwise
	template <class Ind>
	bool readGeneration(std::vector<Ind>& generation, int genNum) {
		std::ifstream file(fileName.c_str());
		if (!file)
			return false;
		
		if (genNum < 0) {
			genNum = analyze<Ind>(-1, -1).generations.back().cardinalNum;
		}
		
		// find genNum generation
		bool found = false;
		int readMode = -1;
		for(;;) {
			if (file.peek() == '#')
				file.ignore(100000, '\n');
			else {
				FileLine<Ind> fline;
				int g = readLine(fline, file, genNum, readMode);
				if (g == ~0) {
					if (found) break;
					else
						throw std::runtime_error("Error in readGeneration - could not read file");
				} else if (g < 0) {
					file.ignore(100000, '\n');
					if (-g == genNum)
						throw std::runtime_error("Error in readGeneration - could not read individual of target generation");
				} else if (g == genNum) {
					found = true;
					generation.push_back(fline.individual);
				} else if (g > genNum) {
					break;
				}
			}
		}
		return found;
	}
	
	// returns false if file not found or not right format, true otherwise
	template <class Ind, class Func>
	bool readGenerations(size_t genNum1, size_t genNum2, FileLine<Ind> fline, Func& func, int readMode = -1) {
		std::cerr << "processing file:\n";
		
		std::ifstream file(fileName.c_str());
		if (!file) {
			std::cerr << "  failed to open " << fileName << "\n";
			return false;
		}
		
		// find genNum generation
		for(;;) {
			if (file.peek() == '#') {
				file.ignore(100000, '\n');
			} else {
				//Ind individual;
				int g = readLine(fline, file, genNum1, readMode);
				
				if (file && ((size_t)g <= genNum2)) {
					std::cerr << "  " << g << "\r";
				} else {
					std::cerr << "  done\n";
					break;
				}
				
				if ((size_t)g >= genNum1)
					func(fline, g);
			}
			if (!file)
				break;
		}
		return true;
	}
	
	template <class Ind>
	void outputGeneration(const std::vector<Ind>& pop, size_t genNum) {
		std::ofstream file;
		if (genNum == 0) {
			file.open(fileName.c_str());
			file <<	"# format of this file:\n" <<
				"# generation_number \t[function_input_vector] \t violation \t[properties] \t[function_output_vector]\n";
		} else {
			file.open(fileName.c_str(), std::ios::app);
		}
		
		for (size_t i = 0; i < pop.size(); ++i) {
			file << genNum << "\t" << pop[i].chromosome << "\t" << pop[i].violation << "\t" 
				<< pop[i].getProperties() << "\t" << pop[i].criteria << "\n";
		}
	}
	
	// function used to log all evaluated individuals
	// function accepts:
	//   the individual
	//   serial number of the evaluation, rank of the processor that did the evaluation 
	//   times of evaluation and total life (from generation until the registration of result)
	template <class Ind>
	void outputEvaluation(const Ind& ind, size_t evalNum, int evalRank, double evalTime, double lifeTime) {
		std::ofstream file;
		if (evalNum == 0) {
			file.open(fileName.c_str());
			file <<	"# format of this file:\n" <<
				"# evaluation_number \t[function_input_vector] \t violation \t[properties] \t[function_output_vector] \t evaluator_rank \t evaluation_time \t life_time \n";
		} else {
			file.open(fileName.c_str(), std::ios::app);
		}
		
		file << evalNum << "\t" << std::setprecision(26) << ind.chromosome;
		file << "\t" << ind.violation << "\t" 
			<< ind.getProperties() << "\t" << ind.criteria << "\t" 
			<< evalRank << "\t" << evalTime << "\t" << lifeTime << "\n";
	}
	
	template <class Ind>
	AnalysisResult analyze(int firstGeneration, int lastGeneration) {
		AnalysisResult result;
		
		std::ifstream file(fileName.c_str());
		std::string line;
		if (!file) {
			result.failedToOpenFile = true;
			result.bigError = true;
			return result;
		}
		
		// analysis
		int readMode = -1;
		while (std::getline(file, line)) {
			// skip empty lines
			if (line.size() == 0)
				continue;
			
			// read comments and header (special type of comment)
			else if (line[0] == '#') {
				std::stringstream strLine(line);
				strLine.ignore(1);
				std::string text;
				strLine >> text;
				// header may specify type of information saved in file
				if (text == "output:") {
					if (readMode == -1)
						readMode = 0;
					
					strLine >> text;
					if (text == "chromosome")
						readMode |= read_mode_chromosome;
					if (text == "violation")
						readMode |= read_mode_violation;
					if (text == "criteria")
						readMode |= read_mode_criteria;
					if (text == "properties")
						readMode |= read_mode_properties;
				}
				continue;
			}
			
			std::istringstream lineStr(line);
			FileLine<Ind> fline;
			bool errors = false;
			int g = readLine(fline, lineStr, firstGeneration, readMode);
			if (g == ~0) {
				result.bigError = true;
				break;
			}
			if (g < 0) {
				g = -g;
				errors = true;
			}
			
			if (result.generations.size() == 0)
				result.generations.push_back(AnalysisResult::Generation(g));
			if (g != result.generations.back().cardinalNum) {
				result.generations.push_back(AnalysisResult::Generation(g));
				std::cerr << "analyzing generation " << g << "\r";
			}
			++result.generations.back().populationSize;
			result.generations.back().errors = errors;
		}		
		result.fileMode = readMode;
		
		for (size_t i = 1; (i < result.generations.size()) && !result.bigError; ++i) {
			result.orderedGenerations &= (result.generations[i].cardinalNum > 
				result.generations[i-1].cardinalNum);
			result.staticPopulationSize &= (result.generations[i].populationSize == 
				result.generations[i-1].populationSize);
			result.errors |= result.generations[i].errors;
		}
		
		if (result.generations.size() == 0) {
			result.bigError = true;
			result.generations.push_back(AnalysisResult::Generation());
		}
		
		std::cerr << "                           \r";
		return result;
	}
};


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// square function
template <class T>
inline T sqr(T a) {
	return a*a;
}


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Pareto domination
template<class T, size_t N>
bool operator< (const Array<T, N>& a, const Array<T, N>& b) {
	bool dom 	= (a[0] < b[0]);
	bool cover 	= (a[0] <= b[0]);
	for (size_t i = 1; i < N; ++i) {
		dom 	|= (a[i] < b[i]);
		cover 	&= (a[i] <= b[i]);
	}
	return dom & cover;
}


template<class T, size_t N>
bool operator<= (const Array<T, N>& a, const Array<T, N>& b) {
	bool ret = (a[0] <= b[0]);
	for (size_t i = 1; i < N; ++i)
		ret &= (a[i] <= b[i]);
	return ret;
}


template<class T>
bool operator< (const std::vector<T>& a, const std::vector<T>& b) {
	assert(a.size() == b.size());
	
	bool dom 	= (a[0] < b[0]);
	bool cover 	= (a[0] <= b[0]);
	for (size_t i = 1; i < a.size(); ++i) {
		dom 	|= (a[i] < b[i]);
		cover 	&= (a[i] <= b[i]);
	}
	return dom & cover;
}


template<class T>
bool operator<= (const std::vector<T>& a, const std::vector<T>& b) {
	assert(a.size() == b.size());
	
	bool ret = (a[0] <= b[0]);
	for (size_t i = 1; i < a.size(); ++i)
		ret &= (a[i] <= b[i]);
	return ret;
}


/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Evklidian distance
template<class Vec>
double vectorDistance(const Vec& a, const Vec& b) {
	assert(a.size() == b.size());
	
	double dist = sqr(a[0] - b[0]);
	for (size_t i = 1; i < a.size(); ++i) 
		dist += sqr(a[i] - b[i]);
		
	return sqrt(dist);
}


/// helper class for IndividualStruct, used for conditional inclusion of Properties
template<class Properties>
struct IncludeProperties {
	Properties properties;
	
	void operator= (const Properties& other) {
		properties = other.properties;
	}
	
	Properties& getProperties() {return properties;}
	const Properties& getProperties() const {return properties;}
	
	static const bool hasProperties = true;
};


template<>
struct IncludeProperties<void> {
	const char* getProperties() const {return "/";}
	
	static const bool hasProperties = false;
};


/// Definition of an individual:
///  - chromosome
///  - criteria
///  - viloation: value > 0 signals a violation of chromosome value and therefore 
/// 	invalid criteria; violation == -1 signals individual is not evaluated yet
///   
///
template<class Chromosome, class Criteria, class Properties = void>
struct IndividualStruc : IncludeProperties<Properties> {
	Chromosome chromosome;
	Criteria criteria;
	// properties is the information about the individual that is not used within the
	// optimization procedure but is returned by the evaluation function and is worth
	// recording for later analysis of the results
	
	double violation;
	
	IndividualStruc() : violation(-1.0) {}
	
	const IndividualStruc& operator= (const IndividualStruc& other) {
		chromosome = other.chromosome;
		criteria = other.criteria;
		violation = other.violation;
		IncludeProperties<Properties>::operator=(other);
		return *this;
	}
	
	bool evaluated() const {
		return violation != -1;
	}
	
	inline friend bool operator< (const IndividualStruc& ind1, const IndividualStruc& ind2) {
		/*
		if (ind1.violation < 0) 
			std::cerr << " !! a < b: a is invalid\n";
		if (ind2.violation < 0) 
			std::cerr << " !! a < b: b is invalid\n";
			*/
		return (ind1.violation < ind2.violation) || 
			((ind1.violation == ind2.violation) && (ind1.criteria < ind2.criteria));
	}
	
	inline friend bool operator<= (const IndividualStruc& ind1, const IndividualStruc& ind2) {
		/*
		if (ind1.violation < 0) 
			std::cerr << " !! a <= b: a is invalid\n";
		if (ind2.violation < 0) 
			std::cerr << " !! a <= b: b is invalid\n";
			*/
		return (ind1.violation < ind2.violation) || 
			((ind1.violation == ind2.violation) && (ind1.criteria <= ind2.criteria));
	}
};

#endif // INDIVIDUAL_H_INCLUDED
