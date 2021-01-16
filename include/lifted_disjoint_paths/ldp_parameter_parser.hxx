#ifndef LDP_PARAMETER_PARSER_HXX
#define LDP_PARAMETER_PARSER_HXX

#include<map>
#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include"ldp_methods.hxx"

namespace LPMP {


class ParametersParser{
  public:
    ParametersParser(){}


    void initFromFile(std::string& fileName,bool searchSolverKeyword){
        std::ifstream data;
        //data.exceptions( std::ifstream::failbit | std::ifstream::badbit );
        try{
            data.open(fileName);
            if (!data){
                throw std::system_error(errno, std::system_category(), "failed to open parameter file "+fileName);
            }
            bool solverPart=false;
            std::string line;
            if(searchSolverKeyword){
                while (!solverPart&&std::getline(data, line) ) {
                    size_t commentPos=line.find_first_of('#');
                    if(commentPos<line.length()){
                        line.erase(commentPos);
                    }
                    if(!line.empty()&&line.find("[SOLVER]")!=line.npos){
                        solverPart=true;
                    }
                }
                if(!solverPart){
                    throw std::runtime_error("Config file does not contain \"[SOLVER]\".  ");
                }
            }
            initFromStream(data);

            data.close();

        }

        catch (std::system_error& er) {
            std::clog << er.what() << " (" << er.code() << ")" << std::endl;
        }

    }
    void initFromFile(std::string& fileName){
        initFromFile(fileName,true);
    }


    template<class STR>
    void initFromStream(STR& data){
        char delim='=';
        std::string line;
        std::vector<std::string> strings;

        bool newSection=false;
        while(!newSection&&std::getline(data, line)){
            size_t commentPos=line.find_first_of('#');
            if(commentPos<line.length()){
                line.erase(commentPos);
            }
            if(line.find("[")==line.npos){
                if(!line.empty()){ //TODO not split all delims, just the first occurence
                    strings=split<>(line,delim);
                    std::string whitespaces (" ");

                    size_t foundLast = strings[0].find_last_not_of(whitespaces);
                    size_t foundFirst=strings[0].find_first_not_of(whitespaces);
                    std::string key=strings[0].substr(foundFirst,foundLast-foundFirst+1);

                    foundLast = strings[1].find_last_not_of(whitespaces);
                    foundFirst=strings[1].find_first_not_of(whitespaces);
                    std::string value=strings[1].substr(foundFirst,foundLast-foundFirst+1);

                    parsedStrings[key]=value;
                }
            }
            else{
                newSection=true;
            }
        }

    }


    std::map<std::string, std::string>& getParsedStrings(){
        return parsedStrings;
    }

    void clearParsedStrings(){
        parsedStrings.clear();
    }

private:
    std::map<std::string,std::string> parsedStrings;

};


}


#endif // LDP_PARAMETER_PARSER_HXX
