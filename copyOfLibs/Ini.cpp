#include "Ini.h"


//#include <limits>
#include <iostream>
#include <fstream>


namespace Ini {

    bool isElementOf(char ch, const std::string& str) {
        for (unsigned int i=0; i<str.length(); i++) {
            if (ch == str[i])
                return true;
        }
        return false;
    }


    Name::Name(const char* name)
        : name_(name) {}


    Name::Name(const std::string& name)
            : name_(name) {}


    bool operator== (const std::string& strName, const Name& name) {
        return strName == name.name_;
    }


    Name::operator const std::string&() const {
        return name_;
    }


    Variable::Variable(const std::string& name, const std::string& val, int lines)
            : name_(name),
            value_(val),
            precedingEmptyLines_(lines) {}


    Section::Section(const std::string& name, int lines)
            : precedingEmptyLines_(lines) {
        // directly add default
        if (name == "default") {
            name_ = name;
        } else {
            name_ = normalizeName(name);
        }
    }


    std::deque<Variable>::const_iterator Section::findName(const Name& name) const {
        for (std::deque<Variable>::const_iterator it = variables_.begin(); it != variables_.end(); ++it) {
            if (it->name_ == name)
                return it;
        }

        return variables_.end();
    }


    std::deque<Variable>::iterator Section::findName(const Name& name) {
        for (std::deque<Variable>::iterator it = variables_.begin(); it != variables_.end(); ++it) {
            if (it->name_ == name)
                return it;
        }

        return variables_.end();
    }


    Name Section::normalizeName(const Name& name) {
        // lose starting whitespace and whitespace between the [ char and the actual section name
        int                        start, count;
        {
            bool                    noBracket = true;
            for (start = 0; start < (int)name.str().length(); ++start) {
                if ( (name.str()[start] == ' ') || (name.str()[start] == '\t') ) {}
                else if (noBracket && (name.str()[start] == '[')) {
                    noBracket = false;
                } else
                    break;
            }
        }

        // lose trailing whitespace and whitespace after the name and before the ] char
        {
            bool                    noBracket = true;
            for (count = (int)name.str().length()-1; count > 0; --count) {
                if ( (name.str()[count] == ' ') || (name.str()[count] == '\t') ) {}
                else if (noBracket && (name.str()[count] == ']')) {
                    noBracket = false;
                } else
                    break;
            }
            count += 1 - start;
        }

        // use normalized name
        return name.str().substr(start, count);
    }


    void Section::compact() {
        // a double loop that leaves only single occurance of each found variable (renames others
        // into "")
        for (int i = 0; i < (int)variables_.size(); ++i) {
            if (variables_[i].name_ == "")
                continue;

            if (variables_[i].name_[0] == ';') {
                variables_[i].name_ = "comment";
                continue;
            }

            for (int j = i+1; j < (int)variables_.size(); ++j) {
                //std::cout << "checking " << i << " <-> " << j << "\n";
                if (variables_[i].name_ == variables_[j].name_) {
                    //std::cout << "cleaning section " << variables_[i].name_ << " " << i << " [" << j << "]\n";
                    variables_[i].name_ = "";
                    break;
                }
            }
        }

        // a loop that cleans up all occurances of variable named ""
        int                        l = 0;
        for (int k = 0; k < (int)variables_.size(); ++k) {
            if (variables_[k].name_ != "") {
                if (k != l)
                    variables_[l] = variables_[k];
                ++l;
            }
        }
        variables_.erase(variables_.begin() + l, variables_.end());
    }


    /// ***********************************************************************************************
	/// class File
	/// 
	File::File(const char* filename)
            : separators_(" \t"),
            filename_(filename),
            state_(fileState_error) {
        std::ifstream              file;
        file.open(filename);
        if (!file.is_open()) {
            state_ = fileState_notFound;
            return;
        } else {
            state_ = fileState_openToRead;
        }

        std::string                line;
        std::string                name, val;
        int                        emptyLines = 0;
        int                        activeSection = 0;
        sections_.push_back(Section("default"));

        while (file.good()) {
            // parse single line of file and create two strings - (variable) name and value
            if (file.eof()) break;
            std::getline(file, line);
            if (line == "") {
                ++emptyLines;
                continue;
            }
            parse(line, name, val);

            // did we just parse a section or a variable?
            if (isSection(name)) {
                // a section
                activeSection = addSection(name, emptyLines);
            } else {
                // a variable
                sections_[activeSection].variables_.push_front(Variable(name, val, emptyLines));
            }
            emptyLines = 0;
            //std::numeric_limits<int>::max()
        }

        if ((sections_.size() == 1) && (sections_[0].size() == 0)) {
            state_ |= fileState_empty;
        }
    }


    bool File::save() const {
        try {
            std::ofstream              file;
            file.open(filename_.c_str());
            if (!file.is_open())
                return false;

            for (unsigned int sec=0; sec<sections_.size(); ++sec) {
                if (sec > 0) {
                    for (int i = 0; i < sections_[sec].precedingEmptyLines_; ++i)
                        file << "\n";
                    file << "[" << sections_[sec].name_ << "]" << std::endl;
                }

                for (std::deque<Variable>::const_reverse_iterator it = sections_[sec].variables_.rbegin();
                        it != sections_[sec].variables_.rend(); ++it) {
                    for (int i = 0; i < it->precedingEmptyLines_; ++i)
                        file << "\n";
                    file << it->name_;
                    if (it->value_ != "")
                        file << " = " << it->value_;
                    file << "\n";
                    if (!file.good())
                        return false;
                }
            }
            file << std::endl;
            if (state_ & fileState_notFound) {
                state_ = fileState_created;
            } else {
                state_ = 0;
            }
        } catch (...) {
            state_ = fileState_error;
            return false;
        }
        return true;
    }


    int File::getSectionNumber(const std::string& name) const {
        for (int sec=1; sec<(int)sections_.size(); ++sec) {
            if (sections_[sec].name_ == name)
                return sec;
        }
        return -1;
    }


    bool File::loadVar(std::string& var, Name name, int section) const {
        if ((section >= (int)sections_.size()) || (section < 0))
            return false;

        std::deque<Variable>::const_iterator   it = sections_[section].findName(name);
        if (it != sections_[section].variables_.end()) {
            var = it->value_;
            return true;
        }
        return false;
    }


    int File::addSection(Name name, int emptyLines) {
        name = Section::normalizeName(name);
        int                        foundSection = findSection(name);
        if (foundSection == -1) {
            // didn't find the section -> create new
            foundSection = (int)sections_.size();
            sections_.push_back(Section(name, emptyLines));
        }
        return foundSection;
    }


    void File::compactSection(int section) {
        if ((section >= (int)sections_.size()) || (section < 0))
            return;

        sections_[section].compact();
    }


    void File::print() {
        for (unsigned int sec=0; sec<sections_.size(); ++sec) {
            std::cout << sections_[sec].name_ << "\n";

            for (std::deque<Variable>::const_iterator it = sections_[sec].variables_.begin();
                    it != sections_[sec].variables_.end(); ++it) {
                std::cout << "   " << it->name_;
                if (it->value_ != "")
                    std::cout << " = " << it->value_;
                std::cout << "\n";
                if (!std::cout.good()) {
                    std::cout.clear();
                    return;
                }
            }
        }
    }


    void File::parse(const std::string& input, std::string& name, std::string& val) const {
        int                        lastChar = 0, firstChar = 0;
        bool                       foundName = false;
        for (unsigned int i=0; i<input.length(); ++i) {
            if ((input[i] == '=') && (!foundName)) {
                name = std::string(input, firstChar, lastChar - firstChar + 1);
                foundName = true;
                firstChar = i+1;
            }

            if (isElementOf(input[i], separators_)) {
                if (i == firstChar)
                    ++firstChar;
            } else {
                lastChar = i;
            }
        }

        if (foundName) {
            val = std::string(input, firstChar, lastChar - firstChar + 1);
        } else {
            name = std::string(input, firstChar, lastChar - firstChar + 1);
            val = "";
        }
    }


    bool File::isSection(const std::string& name) const {
        return ((name[0]=='[') && (name[name.length()-1]==']'));
    }


    int File::findSection(const std::string& name) const {
        for (unsigned int i=0; i<sections_.size(); ++i) {
            if (sections_[i].name_ == name) {
                return i;
            }
        }
        return -1;
    }


    CharSeparator::CharSeparator(const std::string& chars)
            : chars_(chars) {}


    bool CharSeparator::removeFrom(std::istream& inStream) {
        int                        removeCount = 0;
        while (inStream.good()) {
            char                    ch = inStream.peek();
            if (isElementOf(ch, chars_)) {
                inStream.ignore();
                ++removeCount;
            } else
                break;
        }

        return (removeCount > 0);
    }


    bool CharSeparator::removeFirst(std::istream& inStream) {
        return true;
    }

}
