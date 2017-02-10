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
