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

#ifndef ARGUMENTS_H
#define ARGUMENTS_H


#include <string>
#include <sstream>
#include <map>
#include <vector>


class Arguments {
protected:
    int                           num_;
    char**                        args_;

public:
    Arguments(int num, char** args);

    const char* programName() const;
    int size() const;
    const char* argument(int num) const;
    const char* operator[] (int num) const;
};


// ************************************************************************************************
// class ArgumentsCStyle
/** 
	Colon Style arguments (e.g. -in:input.txt -out:screen)
**/
class ArgumentsCStyle : public Arguments {
    typedef std::map<std::string, std::string> StrStrMap;
    typedef std::map<std::string, std::string>::const_iterator StrStrMapIt;
    StrStrMap                     argMap_;

public:
    ArgumentsCStyle(int num, char** args, char colon=':');

    const std::string& operator[] (const std::string& name) const;

    template<class T>
    bool setVar(T& var, const std::string& name) {
        StrStrMapIt                it = argMap_.find(name);
        if (it != argMap_.end()) {
            std::stringstream       stream(it->second);
            stream >> var;
            return true;
        } else
            return false;
    }

    bool setVar(std::string& var, const std::string& name) {
        StrStrMapIt                it = argMap_.find(name);
        if (it != argMap_.end()) {
            var = it->second;
            return true;
        } else
            return false;
    }
};


/// helper function for DashArguments
template<class Var>
bool extractVar(const std::string& inputString, Var& var) {
	std::stringstream stream(inputString);
	stream >> var;
	return bool(stream);
}

/// overload of basic extractVar function
template<class Ch, class ChTr>
bool extractVar(const std::string& inputString, std::basic_string<Ch, ChTr>& var) {
	var = inputString;
	return true;
}


// ************************************************************************************************
// class DashArguments
/**
    Class for parsing arguments of the form: -var1 value1 -var2 value2.1 value2.2 ... freeArguments
    Features:
      - Variable is formed of a dash argument (var name) and non dash argument (value);
      - It allows several vars to have identical names (multiple values per var)
      - It allows for flags (variables with no value)
      - FreeArguments are arguments without the dash, that are not values (that do not directly follow a dash argument)
	@todo add option to make arguments with equality sign
	@warning when checking for individual arguments don't forget to write dash as the part of the argument name
	@code DashArguments.setVar(a, "-a");
**/
class DashArguments : public Arguments {
	typedef std::multimap<std::string, std::string> StrStrMap;
	typedef std::multimap<std::string, std::string>::iterator StrStrMapIt;
	typedef std::multimap<std::string, std::string>::const_iterator constStrStrMapIt;
	StrStrMap					dashArgs;
	std::vector<std::string>	freeArgs;

public:
    /// in ctor, the whole ?.ini file is read, parsed and (variable, value) tuples are stored in memory
	DashArguments(int num, char** args);

	const std::string& operator[] (const std::string& name) const;

    /// This is the function the whole class is built around. It reads variable value from memory
    /// returns false if the variable cannot be found, true otherwise (even if it has no value)
	template<class T>
	bool setVar(T& var, const std::string& name) {
		constStrStrMapIt it = dashArgs.find(name);
		if (it != dashArgs.end()) {
			return extractVar(it->second, var);
		} else
			return false;
	}

	/// Dash arguments allows multiple arguments with same name (and different values)
	/// This function allows exploring all the values of a single argument name
	/// example in program call: myprog.exe -out text -out graph
	/// example in function call: 
	///	  - @code vector<string> outputs;
	///	  - @code setVars("-out", outputs)
	template<class T>
	bool setVars(std::vector<T>& vars, const std::string& name) {
		std::pair<StrStrMapIt, StrStrMapIt> itPair = dashArgs.equal_range(name);
		bool extractionOk = true;
		for (StrStrMapIt it = itPair.first; it != itPair.second; ++it) {
			T var;
			bool ok = extractVar(it->second, var);
			if (ok)
				vars.push_back(var);
			extractionOk = extractionOk && ok;
		}
		return extractionOk && (itPair.first != itPair.second);
	}
	
	/// check if flag is set (variable exists but may have no value)
	bool isSet(const std::string& name) {
		constStrStrMapIt it = dashArgs.find(name);
		return (it != dashArgs.end());
	}
	
	/// returns arguments that were not parsed (all of those that do not start with a dash)
	const std::vector<std::string>& freeArguments() const {
		return freeArgs;
	}
};


#endif //ARGUMENTS_H
