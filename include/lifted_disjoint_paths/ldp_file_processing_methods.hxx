#ifndef LDP_FILE_PROCESSING_METHODS_HXX
#define LDP_FILE_PROCESSING_METHODS_HXX

#include<string>
#include<vector>
#include<iostream>
#include<fstream>


namespace LPMP {


template<class T=char>
std::vector<std::string> split(
        std::string inputString, T delim) {
    size_t occurence = 0;
    size_t newOccurence = 0;
    std::vector<std::string> strings;
    while (newOccurence < inputString.size()) {
        newOccurence = std::min(inputString.find_first_of(delim, occurence),
                                inputString.size());

        std::string newString(inputString, occurence, newOccurence - occurence);
        strings.push_back(newString);
        newOccurence = newOccurence + 1;
        occurence = newOccurence;
    }

    return strings;
}




template<class T=size_t>
std::vector<std::vector<T>> readLines(std::string inputFileName, char delim) {
    std::ifstream inputFile;
    std::vector<std::vector<T>> outputs;
    try{

        inputFile.open(inputFileName);
        if (!inputFile){
            throw std::system_error(errno, std::system_category(), "failed to open file with interval solution "+inputFileName);
        }

        std::string line;
        size_t lineCounter=0;
        std::getline(inputFile, line);
        std::vector<std::string> strings;

        while (std::getline(inputFile, line) && !line.empty()) {
            strings=split(line,delim);
            size_t length=strings.size();
            std::vector<T> parsedLine(length);
            for (int i = 0; i < length; ++i) {
                parsedLine[i]=std::stoul(strings[i]);
            }
            outputs.push_back(parsedLine);

        }

        inputFile.close();


    }

    catch (std::system_error& er) {
        std::clog << er.what() << " (" << er.code() << ")" << std::endl;

    }
    return outputs;


}



inline void writeOutputToFile(const std::vector<std::vector<size_t>>& paths,std::string outputFileName){


    std::ofstream file;
    file.open(outputFileName);


    for (int i = 0; i < paths.size(); ++i) {
        //std::cout<<"output path "<<i<<std::endl;
        for (int j = 0; j < paths[i].size(); ++j) {
            size_t v=paths[i][j];
            file<<v<<" ";
            //	labels[v]=i+1;
        }
        file<<std::endl;
    }


    file.close();

}



}

#endif // LDP_FILE_PROCESSING_METHODS_HXX
