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

#include <string>
#include <deque>
#include <sstream>
#include <iostream>
#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include <functional>
#include <utility>
#endif //__GXX_EXPERIMENTAL_CXX0X__


#ifndef INI_H
#define INI_H
namespace Ini {

    bool isElementOf(char ch, const std::string& str);


    class Name {
        std::string                name_;

    public:
        Name(const char* name);
        Name(const std::string& name);
        friend bool operator== (const std::string& strName, const Name& name);
        operator const std::string&() const;

        inline std::string& str() {
            return name_;
        }

        inline const std::string& str() const {
            return name_;
        }
    };


    struct Variable {
        std::string                name_;
        std::string                value_;
        int                        precedingEmptyLines_;

        Variable(const std::string& name, const std::string& val, int lines = 0);
    };


    struct Section {
        std::string                name_;
        std::deque<Variable>       variables_;
        int                        precedingEmptyLines_;

	public:
        Section(const std::string& name, int lines = 0);
        std::deque<Variable>::const_iterator findName(const Name& name) const;
        std::deque<Variable>::iterator findName(const Name& name);
        static Name normalizeName(const Name& name);
        int size() const  {
            return (int)variables_.size();
        }

        void compact();
    };


    class CharSeparator {
        std::string                chars_;
    public:
        CharSeparator(const std::string& chars = " ;,\t");
        bool removeFrom(std::istream& inStream);
        bool removeFirst(std::istream& inStream);
    };


	/// ***********************************************************************************************
	/// class File
	/// 
	/// The only class of Ini, that user should work with directly (creating variables of)
	/// 
    class File {
        std::deque<Section>        sections_;
        std::string                filename_;
		
        std::string                separators_;
		mutable int						state_;
		
    public:
        enum FileState {
            fileState_error = 1,
            fileState_openToWrite = 2,
            fileState_openToRead = 4,
            fileState_empty = 8,
            fileState_created = 16,			// file did not exist before but it has been created
            fileState_notFound = 32
        };
		
    public:
        File(const char* filename);
		
        // was file not found?
        bool notFound() const {
            return ((state_ & fileState_notFound) > 0);
        }
		
        // is file empty?
        bool empty() const {
            return ((state_ & fileState_empty) > 0);
        }
		
        // ok when there is no error (file not found is not reported here)
        bool ok() const {
            return ((state_ & fileState_error) == 0);
        }
        
        operator bool() const {
        	return ok() && !notFound();
        }
		
        bool save() const;
		
        // section handlers
        // default section has number 0, result -1 marks an error
        int getSectionNumber(const std::string& name) const;
		
        // variable handlers
        template<class T>
        bool loadVar(T& var, Name name, int section = 0) const {
            if ((section >= (int)sections_.size()) || (section < 0))
                return false;
			
            std::deque<Variable>::const_iterator    it = sections_[section].findName(name);
            if (it != sections_[section].variables_.end()) {
                std::istringstream tempSS(it->value_);
                tempSS >> var;
                return true;
            }
            return false;
        }
		
        bool loadVar(std::string& var, Name name, int section = 0) const;

#ifdef __GXX_EXPERIMENTAL_CXX0X__
        // note: it seems necessary to call LoadArray with eplicitly set template parameter T,
        // e.g. loadArray<double>([myArray&](double& d, int){myArray.push_back(d);});
        template<class T>
        bool loadArray(std::function<void(T, int)> f, Name name, const std::string separators = " ;,\t", int section = 0) const {
            if ((section >= (int)sections_.size()) || (section < 0))
                return false;
			
            std::deque<Variable>::const_iterator    it = sections_[section].findName(name);
                
            if (it != sections_[section].variables_.end()) {
                CharSeparator separator(separators); // TODO : other separator types
                std::istringstream inStream(it->value_);
                
                if (!separator.removeFirst(inStream))  // does not do anything for CharSeparator...
                    return false;

                for (int cnt=0; inStream; ++cnt) {
                    typename std::remove_const<typename std::remove_reference<T>::type>::type temp;
                    inStream >> temp;
                    f(temp, cnt);
                    if (!separator.removeFrom(inStream))
                        break;
                }
            }
            return false;
        }
#endif //__GXX_EXPERIMENTAL_CXX0X__ 
		
        int addSection(Name name, int emptyLines = 0);
        void compactSection(int section);
        template<class T>
        void storeVar(const T& var, Name name, int section = 0) {
            if (section >= sections_.size())
                return;
			
            std::stringstream                   tempStream;
            tempStream << var;
            std::deque<Variable>::iterator      it = sections_[section].findName(name);
            if (it != sections_[section].variables_.end()) {
                it->value_ = tempStream.str();
            } else {
                sections_[section].variables_.push_front(Variable(name, tempStream.str()));
            }
        }
		
        void print();
		
    protected:
        void parse(const std::string& input, std::string& name, std::string& val) const;
        bool isSection(const std::string& name) const;
        int findSection(const std::string& name) const;
    };


    template <class T, class Container>
    struct addToStdContainer {
        void operator() (Container& cont, const T& element) {
            cont.push_back(element);
        }
    };


    template <class T, class ArrayT, class AddFunc = addToStdContainer<T, ArrayT>, class Separator = CharSeparator>
    class ArrayReader {
        ArrayT&                    target_;
        unsigned int               limit_;
        Separator                  separator_;
        AddFunc							addFunc_;

    public:
        ArrayReader(ArrayT& target, unsigned int limit = -1)
                : target_(target),
                limit_(limit) {}

        friend std::istream& operator>>(std::istream& inStream, ArrayReader& reader) {
            if (!reader.separator_.removeFirst(inStream))
                return inStream;

            while (inStream) {
                T                    temp;
                inStream >> temp;
                reader.addFunc_(reader.target_, temp);
                if (!reader.separator_.removeFrom(inStream))
                    break;

                if (reader.limit_ != 0) {
                    --reader.limit_;
                    if (reader.limit_ == 0)
                        break;
                }
            }

            return inStream;
        }
    };

}
#endif // INI_H
