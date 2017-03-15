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

#include "Arguments.h"


Arguments::Arguments(int num, char** args)
		: num_(num - 1),
        args_(args) 
{
}


const char* Arguments::programName() const {
    return args_[0];
}


int Arguments::size() const {
    return num_;
}


const char* Arguments::argument(int num) const {
    return args_[num + 1];
}


const char* Arguments::operator[] (int num) const {
    return args_[num + 1];
}


ArgumentsCStyle::ArgumentsCStyle(int num, char** args, char colon)
        : Arguments(num, args) 
{
    for (int i=0; i<num_; ++i) {
    	if (Arguments::operator[](i) == 0)
			continue;
        std::string                temp = Arguments::operator[](i);
        for (unsigned int j=0; j<temp.size(); ++j) {
            if (temp[j] == colon) {
                argMap_[std::string(temp, 0, j)] = std::string(temp, j+1);
                break;
            }
        }
    }
}


const std::string& ArgumentsCStyle::operator[] (const std::string& name) const {
    static std::string            notFound("");
    StrStrMapIt                   it = argMap_.find(name);
    if (it != argMap_.end())
        return it->second;
    else
        return notFound;
}


DashArguments::DashArguments(int num, char** args) 
   : Arguments(num, args)
{
	StrStrMapIt lastDashArg = dashArgs.end();
	for (int i=0; i<num_; ++i) {
		std::string temp = Arguments::operator[](i);
	    if ((temp.size() > 1) && (temp[0] == '-') && (!isdigit(temp[1]))) {
	    	lastDashArg = dashArgs.insert(std::make_pair(temp, std::string()));
	    } else {
			if ((dashArgs.size() > 0) && (lastDashArg->second == "")) {
				lastDashArg->second = temp;
			} else {
				freeArgs.push_back(temp);
			}
	    }
	}
}


const std::string& DashArguments::operator[] (const std::string& name) const {
	static std::string            notFound("");
	constStrStrMapIt              it = dashArgs.find(name);
	if (it != dashArgs.end())
		return it->second;
	else
		return notFound;
}
