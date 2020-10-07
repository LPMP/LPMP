/*
 * inputConfix.hxx
 *
 *  Created on: Jun 27, 2019
 *      Author: fuksova
 */

#ifndef INCLUDE_DISJOINT_PATHS_INPUTCONFIG_HXX_
#define INCLUDE_DISJOINT_PATHS_INPUTCONFIG_HXX_

#include <string>
//#include <sstream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <limits>
#include "disjoint-paths/parametersParser.hxx"

namespace LPMP{
namespace lifted_disjoint_paths {

template<class T=size_t>
std::vector<std::string> split(
		std::string inputString, char delim) {
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


template<class T = size_t>
class ConfigDisjoint {
public:

	~ConfigDisjoint(){
		std::cout<<"config disjoint destructor"<<std::endl;
		infoFile()<<"config disjoint destructor"<<std::endl;
		(*pInfoFile).close();
		delete pInfoFile;
	}


    ConfigDisjoint(const std::string& inputFileName);
    ConfigDisjoint(std::map<std::string,std::string>& parameters);


	const std::string& getGraphFileName() const {
		return graphFileName;
	}

	size_t getMaxTimeBase() const {
		return maxTimeBase;
	}

	size_t getMaxTimeFrame() const {
		return maxTimeFrame;
	}

	const size_t getMaxTimeLifted() const {
		//std::cout<<"getting max time lifted "<<maxTimeLifted<<std::endl;
		return maxTimeLifted;
	}



	const std::string& getTimeFileName() const {
		return timeFileName;
	}


	double getBaseUpperThreshold() const {
		return baseUpperThreshold;
	}

	size_t getDenseTimeLifted() const {
		return denseTimeLifted;
	}

	size_t getKnnK() const {
		return knnK;
	}

	size_t getKnnTimeGap() const {
		return knnTimeGap;
	}

	size_t getLongerIntervalLifted() const {
		return longerIntervalLifted;
	}

	double getNegativeThresholdLifted() const {
		return negativeThresholdLifted;
	}

	double getPositiveThresholdLifted() const {
		return positiveThresholdLifted;
	}


	bool isRestrictFrames() const {
		return restrictFrames;
	}

	double getInputCost() const {
		return inputCost;
	}

	double getOutputCost() const {
		return outputCost;
	}

	bool isSparsify() const {
		return sparsify;
	}

	bool isDebugOutputFiles() const {
		return debugOutputFiles;
	}



//	void setMaxTimeLifted(size_t maxTimeLifted) {
//		this->maxTimeLifted = maxTimeLifted;
//	}


	size_t getMaxTimeGapComplete() const {
		return maxTimeGapComplete;
	}


	std::ofstream& infoFile() const {
		std::ofstream & infoF=*pInfoFile;
		return infoF;
	}

    std::stringstream& getControlOutput(){
          return controlOutput;
      }


    void writeControlOutput(){
           std::cout<<controlOutput.str();

//            if(controlOutputFiles){
//                infoFile()<<controlOutput.str();
//                infoFile().flush();
//            }
            controlOutput=std::stringstream();
        }


	//ConfigDisjoint<>& operator=(const ConfigDisjoint<>&);


private:
	ConfigDisjoint<T>(const ConfigDisjoint<T>&);
	ConfigDisjoint();

    void init(std::map<std::string,std::string>& parameters);

    //std::pair<std::string,std::string> parseLine(std::string line,char delim);
	//std::ofstream * fileForGeneralOutputs;
    mutable std::ofstream* pInfoFile;
    std::string inputFileName;  //file with 3 paths: graph, time information, parameters
    std::string outputFileName;   //file (prefix) where outputs are going to be stored

	std::string graphFileName;  //file with input graph structure
	std::string timeFileName;  //file with input time frames information

	std::string paramsFileName; //file where input parameters are saved for future reference
	std::string inputParamsFileName; //file where input parameters are read
	bool debugOutputFiles;  //use debugging output files

	double inputCost;   //cost of input edges (starting a track)
	double outputCost;  //cost of output edges (terminating a track)

	bool sparsify;     //use graph sparsification?
	bool restrictFrames;  //is maximal number of frames given
	size_t maxTimeFrame;   //maximal time frame to be used
	//size_t smallIntervals;


	size_t maxTimeBase;  //max timegap for base edges
	double baseUpperThreshold;  //upper cost thereshold for base edges
	size_t knnTimeGap;   //time gap for base edges to be added densely (all K best to every node in every layer)
	size_t knnK;         //Parameter K for knn

	//bool automaticLifted;
	double negativeThresholdLifted;  //negative cost threshold for lifted edges (lifted edges with cost close to zero are not included)
	double positiveThresholdLifted;  //positive cost threshold for lifted edges (lifted edges with cost close to zero are not included)
	size_t maxTimeLifted;            //max time gap for lifted edges
	size_t denseTimeLifted;          //timegap for lifted edges to be added densely
	size_t longerIntervalLifted;    //time frame skip for adding lifted edges sparsely

	size_t maxTimeGapComplete;     //max time gap of edges to read

     std::stringstream controlOutput;


};

//template<class T>
//inline std::pair<std::string,std::string> ConfigDisjoint<T>::parseLine(std::string line,char delim){
//	std::vector<std::string> strings;
//	strings=split<>(line,delim);
//	std::string whitespaces (" ");

//	size_t foundLast = strings[0].find_last_not_of(whitespaces);
//	size_t foundFirst=strings[0].find_first_not_of(whitespaces);
//	std::string key=strings[0].substr(foundFirst,foundLast-foundFirst+1);

//	foundLast = strings[1].find_last_not_of(whitespaces);
//	foundFirst=strings[1].find_first_not_of(whitespaces);
//	std::string value=strings[1].substr(foundFirst,foundLast-foundFirst+1);

//	return std::pair<std::string,std::string>(key,value);
//}

template<class T>
inline ConfigDisjoint<T>::ConfigDisjoint(const std::string &inputFileName){
    disjointPaths::ParametersParser parser;
    std::string fileName=inputFileName;
    parser.initFromFile(fileName,false);
    std::map<std::string,std::string>& pathsParameters=parser.getParsedStrings();
    if(pathsParameters.count("INPUT_PARAMS")>0){
        std::string paramsFileName=pathsParameters["INPUT_PARAMS"];
        parser.initFromFile(paramsFileName,true);
    }
    init(pathsParameters);


}

template<class T>
inline ConfigDisjoint<T>::ConfigDisjoint(std::map<std::string,std::string>& parameters){
    init(parameters);
}

template<class T>
inline void ConfigDisjoint<T>::init(std::map<std::string,std::string>& parameters){


    if(parameters.count("INPUT_FRAMES")>0){
        timeFileName=parameters["INPUT_FRAMES"];
    }
    else{
        timeFileName="problemDesc_frames";
    }
    std::cout<<"input frames "<<timeFileName<<std::endl;
    //paramsFile<<"input frames "<<timeFileName<<std::endl;


    if(parameters.count("INPUT_GRAPH")>0){
        graphFileName=parameters["INPUT_GRAPH"];
    }
    else{
        throw std::runtime_error("File with the input graph was not specified in the config file");
    }
    std::cout<<"input graph "<<graphFileName<<std::endl;


//    if(pathParameters.count("INPUT_PARAMS")>0){
//        inputParamsFileName=pathParameters["INPUT_PARAMS"];
//    }
//    else{
//        inputParamsFileName="";
//    }
//    std::cout<<"input parameters file "<<inputParamsFileName<<std::endl;

    if(parameters.count("OUTPUT_PATH")==0){
        parameters["OUTPUT_PATH"]=graphFileName.substr(0,graphFileName.find_last_of('/'));
    }

    if(parameters.count("OUTPUT_PREFIX")==0){
        parameters["OUTPUT_PREFIX"]=parameters["OUTPUT_PATH"]+"_output";
    }


    outputFileName=parameters["OUTPUT_PATH"]+parameters["OUTPUT_PREFIX"];
    paramsFileName=outputFileName+"-params.txt";
    std::ofstream paramsFile(paramsFileName.data(),std::ofstream::out);


    paramsFile<<"input graph "<<graphFileName<<std::endl;
    paramsFile<<"output files "<<outputFileName<<std::endl;

    std::string infoFileName=outputFileName+"-info.txt";
    pInfoFile=new std::ofstream(infoFileName.data(),std::ofstream::out);

    std::cout<<"output files "<<outputFileName<<std::endl;


    if(parameters.count("DEBUG_OUTPUT_FILES")>0){
        debugOutputFiles=std::stoi(parameters["DEBUG_OUTPUT_FILES"]);
    }
    else{
        debugOutputFiles=0;
    }
    std::cout<<"debug output files "<<debugOutputFiles<<std::endl;
    paramsFile<<"debug output files "<<debugOutputFiles<<std::endl;


    if(parameters.count("SPARSIFY")>0){
        sparsify=std::stoi(parameters["SPARSIFY"]);
    }
    else{
        sparsify=1;
    }
    std::cout<<"sparsify "<<sparsify<<std::endl;
    paramsFile<<"sparsify "<<sparsify<<std::endl;


    if(parameters.count("MAX_TIMEGAP")>0){
        maxTimeFrame=std::stoul(parameters["MAX_TIMEGAP"]);
        restrictFrames=1;
    }
    else{
        restrictFrames=0;
        maxTimeFrame=std::numeric_limits<size_t>::max();
    }
    std::cout<<"max time frame "<<maxTimeFrame<<std::endl;
    paramsFile<<"max time frame "<<maxTimeFrame<<std::endl;

    if(parameters.count("INPUT_COST")>0){
        inputCost=std::stod(parameters["INPUT_COST"]);
    }
    else{
        inputCost=0;
    }
    std::cout<<"input cost "<<inputCost<<std::endl;
    paramsFile<<"input cost "<<inputCost<<std::endl;

    if(parameters.count("OUTPUT_COST")>0){
        outputCost=std::stod(parameters["OUTPUT_COST"]);
    }
    else{
        outputCost=0;
    }
    std::cout<<"output cost "<<outputCost<<std::endl;
    paramsFile<<"output cost "<<outputCost<<std::endl;

    if(parameters.count("MAX_TIMEGAP_BASE")>0){
        maxTimeBase=std::stoul(parameters["MAX_TIMEGAP_BASE"]);
    }
    else{
        maxTimeBase=60;
    }
    std::cout<<"max time gap base "<<maxTimeBase<<std::endl;
    paramsFile<<"max time gap base "<<maxTimeBase<<std::endl;

    if(parameters.count("KNN_GAP")>0){
        knnTimeGap=std::stoul(parameters["KNN_GAP"]);
    }
    else{
        knnTimeGap=3;
    }
    std::cout<<"KNN time gap "<<knnTimeGap<<std::endl;
    paramsFile<<"KNN time gap "<<knnTimeGap<<std::endl;

    if(parameters.count("KNN_K")>0){
        knnK=std::stoul(parameters["KNN_K"]);
    }
    else{
        knnK=3;
    }
    std::cout<<"KNN k "<<knnK<<std::endl;
    paramsFile<<"KNN k "<<knnK<<std::endl;

    if(parameters.count("BASE_THRESHOLD")>0){
        baseUpperThreshold=std::stod(parameters["BASE_THRESHOLD"]);
    }
    else{
        baseUpperThreshold=0;
    }
    std::cout<<"base upper threshold "<<baseUpperThreshold<<std::endl;
    paramsFile<<"base upper threshold "<<baseUpperThreshold<<std::endl;


    if(parameters.count("MAX_TIMEGAP_LIFTED")>0){
        maxTimeLifted=std::stoul(parameters["MAX_TIMEGAP_LIFTED"]);
    }
    else{
        maxTimeLifted=60;
    }
    std::cout<<"max time gap lifted "<<maxTimeLifted<<std::endl;
    paramsFile<<"max time gap lifted "<<maxTimeLifted<<std::endl;

    if(parameters.count("DENSE_TIMEGAP_LIFTED")>0){
        denseTimeLifted=std::stoul(parameters["DENSE_TIMEGAP_LIFTED"]);
    }
    else{
        denseTimeLifted=20;
    }
    std::cout<<"dense time gap lifted "<<denseTimeLifted<<std::endl;
    paramsFile<<"dense time gap lifted "<<denseTimeLifted<<std::endl;

    if(parameters.count("NEGATIVE_THRESHOLD_LIFTED")>0){
        negativeThresholdLifted=std::stod(parameters["NEGATIVE_THRESHOLD_LIFTED"]);
    }
    else{
        negativeThresholdLifted=-1;
    }
    std::cout<<"negative threshold lifted "<<negativeThresholdLifted<<std::endl;
    paramsFile<<"negative threshold lifted "<<negativeThresholdLifted<<std::endl;

    if(parameters.count("POSITIVE_THRESHOLD_LIFTED")>0){
        positiveThresholdLifted=std::stod(parameters["POSITIVE_THRESHOLD_LIFTED"]);
    }
    else{
        positiveThresholdLifted=1;
    }
    std::cout<<"positive threshold lifted "<<positiveThresholdLifted<<std::endl;
    paramsFile<<"positive threshold lifted "<<positiveThresholdLifted<<std::endl;

    if(parameters.count("LONGER_LIFTED_INTERVAL")>0){
        longerIntervalLifted=std::stoul(parameters["LONGER_LIFTED_INTERVAL"]);
    }
    else{
        longerIntervalLifted=4;
    }
    std::cout<<"longer interval lifted "<<longerIntervalLifted<<std::endl;
    paramsFile<<"longer interval lifted "<<longerIntervalLifted<<std::endl;


    if(parameters.count("MAX_TIMEGAP_COMPLETE")>0){
        maxTimeGapComplete=std::stoul(parameters["MAX_TIMEGAP_COMPLETE"]);
    }
    else{
        maxTimeGapComplete=maxTimeFrame;
    }
    std::cout<<"max time gap complete "<<maxTimeGapComplete<<std::endl;
    paramsFile<<"max time gap complete "<<maxTimeGapComplete<<std::endl;


    paramsFile.close();




}

//template<class T>
//inline ConfigDisjoint<T>::ConfigDisjoint(const std::string &inputFileName){

//    std::cout<<"file name "<<inputFileName<<std::endl;
//    char delim='=';
// 	std::ifstream pathsData(inputFileName);
//    std::cout<<"opened"<<std::endl;
//	std::string line;
//	std::vector<std::string> strings;
//	std::map<std::string,std::string> pathParameters;

//	while(std::getline(pathsData, line) ){
//		size_t commentPos=line.find_first_of('#');
//		if(commentPos<line.length()){
//			line.erase(commentPos);
//		}
//		if(!line.empty()){ //TODO not split all delims, just the first occurence
//			std::pair<std::string,std::string> parsed=parseLine(line,delim);
//			pathParameters[parsed.first]=parsed.second;
//		}
//	}


//	if(pathParameters.count("INPUT_FRAMES")>0){
//		timeFileName=pathParameters["INPUT_FRAMES"];
//	}
//	else{
//		timeFileName="problemDesc_frames";
//	}
//	std::cout<<"input frames "<<timeFileName<<std::endl;
//	//paramsFile<<"input frames "<<timeFileName<<std::endl;


//	if(pathParameters.count("INPUT_GRAPH")>0){
//		graphFileName=pathParameters["INPUT_GRAPH"];
//	}
//	else{
//		throw std::runtime_error("File with the input graph was not specified in the config file");
//	}
//	std::cout<<"input graph "<<graphFileName<<std::endl;


//	if(pathParameters.count("INPUT_PARAMS")>0){
//		inputParamsFileName=pathParameters["INPUT_PARAMS"];
//	}
//	else{
//		inputParamsFileName="";
//	}
//	std::cout<<"input parameters file "<<inputParamsFileName<<std::endl;

//	if(pathParameters.count("OUTPUT_PATH")==0){
//		pathParameters["OUTPUT_PATH"]=graphFileName.substr(0,graphFileName.find_last_of('/'));
//	}

//	if(pathParameters.count("OUTPUT_PREFIX")==0){
//		pathParameters["OUTPUT_PREFIX"]=pathParameters["OUTPUT_PATH"]+"_output";
//	}


//	outputFileName=pathParameters["OUTPUT_PATH"]+pathParameters["OUTPUT_PREFIX"];
//	paramsFileName=outputFileName+"-params.txt";
//	std::ofstream paramsFile(paramsFileName.data(),std::ofstream::out);


//	paramsFile<<"input graph "<<graphFileName<<std::endl;
//	paramsFile<<"output files "<<outputFileName<<std::endl;

//	std::string infoFileName=outputFileName+"-info.txt";
//	pInfoFile=new std::ofstream(infoFileName.data(),std::ofstream::out);

//	std::cout<<"output files "<<outputFileName<<std::endl;


//	std::map<std::string,std::string> parameters;

//	if(!inputParamsFileName.empty()){
//		std::ifstream data(inputParamsFileName);
//		size_t lineCounter=0;
//		std::vector<size_t> currentGroup;
//		bool solverPart=false;
//		while (!solverPart&&std::getline(data, line) ) {
//			size_t commentPos=line.find_first_of('#');
//			if(commentPos<line.length()){
//				line.erase(commentPos);
//			}
//			if(!line.empty()&&line.find("[SOLVER]")!=line.npos){
//				solverPart=true;
//			}
//		}
//		if(!solverPart){
//			throw std::runtime_error("Config file does not contain \"[SOLVER]\".  ");
//		}
//		bool newSection=false;
//		while(!newSection&&std::getline(data, line)){
//			size_t commentPos=line.find_first_of('#');
//			if(commentPos<line.length()){
//				line.erase(commentPos);
//			}
//			if(line.find("[")==line.npos){
//				if(!line.empty()){ //TODO not split all delims, just the first occurence
//					std::pair<std::string,std::string> parsed=parseLine(line,delim);
//					parameters[parsed.first]=parsed.second;
//				}
//			}
//			else{
//				newSection=true;
//			}
//		}

//		for(auto itMap=parameters.begin();itMap!=parameters.end();itMap++){
//			std::cout<<itMap->first<<"-"<<itMap->second<<"-"<<std::endl;
//		}

//	}

//	if(parameters.count("DEBUG_OUTPUT_FILES")>0){
//		debugOutputFiles=std::stoi(parameters["DEBUG_OUTPUT_FILES"]);
//	}
//	else{
//		debugOutputFiles=0;
//	}
//	std::cout<<"debug output files "<<debugOutputFiles<<std::endl;
//	paramsFile<<"debug output files "<<debugOutputFiles<<std::endl;


//	if(parameters.count("SPARSIFY")>0){
//		sparsify=std::stoi(parameters["SPARSIFY"]);
//	}
//	else{
//		sparsify=1;
//	}
//	std::cout<<"sparsify "<<sparsify<<std::endl;
//	paramsFile<<"sparsify "<<sparsify<<std::endl;


//	if(parameters.count("MAX_TIMEGAP")>0){
//		maxTimeFrame=std::stoul(parameters["MAX_TIMEGAP"]);
//		restrictFrames=1;
//	}
//	else{
//		restrictFrames=0;
//		maxTimeFrame=std::numeric_limits<size_t>::max();
//	}
//	std::cout<<"max time frame "<<maxTimeFrame<<std::endl;
//	paramsFile<<"max time frame "<<maxTimeFrame<<std::endl;

//	if(parameters.count("INPUT_COST")>0){
//		inputCost=std::stod(parameters["INPUT_COST"]);
//	}
//	else{
//		inputCost=0;
//	}
//	std::cout<<"input cost "<<inputCost<<std::endl;
//	paramsFile<<"input cost "<<inputCost<<std::endl;

//	if(parameters.count("OUTPUT_COST")>0){
//		outputCost=std::stod(parameters["OUTPUT_COST"]);
//	}
//	else{
//		outputCost=0;
//	}
//	std::cout<<"output cost "<<outputCost<<std::endl;
//	paramsFile<<"output cost "<<outputCost<<std::endl;

//	if(parameters.count("MAX_TIMEGAP_BASE")>0){
//		maxTimeBase=std::stoul(parameters["MAX_TIMEGAP_BASE"]);
//	}
//	else{
//		maxTimeBase=60;
//	}
//	std::cout<<"max time gap base "<<maxTimeBase<<std::endl;
//	paramsFile<<"max time gap base "<<maxTimeBase<<std::endl;

//	if(parameters.count("KNN_GAP")>0){
//		knnTimeGap=std::stoul(parameters["KNN_GAP"]);
//	}
//	else{
//		knnTimeGap=3;
//	}
//	std::cout<<"KNN time gap "<<knnTimeGap<<std::endl;
//	paramsFile<<"KNN time gap "<<knnTimeGap<<std::endl;

//	if(parameters.count("KNN_K")>0){
//		knnK=std::stoul(parameters["KNN_K"]);
//	}
//	else{
//		knnK=3;
//	}
//	std::cout<<"KNN k "<<knnK<<std::endl;
//	paramsFile<<"KNN k "<<knnK<<std::endl;

//	if(parameters.count("BASE_THRESHOLD")>0){
//		baseUpperThreshold=std::stod(parameters["BASE_THRESHOLD"]);
//	}
//	else{
//		baseUpperThreshold=0;
//	}
//	std::cout<<"base upper threshold "<<baseUpperThreshold<<std::endl;
//	paramsFile<<"base upper threshold "<<baseUpperThreshold<<std::endl;


//	if(parameters.count("MAX_TIMEGAP_LIFTED")>0){
//		maxTimeLifted=std::stoul(parameters["MAX_TIMEGAP_LIFTED"]);
//	}
//	else{
//		maxTimeLifted=60;
//	}
//	std::cout<<"max time gap lifted "<<maxTimeLifted<<std::endl;
//	paramsFile<<"max time gap lifted "<<maxTimeLifted<<std::endl;

//	if(parameters.count("DENSE_TIMEGAP_LIFTED")>0){
//		denseTimeLifted=std::stoul(parameters["DENSE_TIMEGAP_LIFTED"]);
//	}
//	else{
//		denseTimeLifted=20;
//	}
//	std::cout<<"dense time gap lifted "<<denseTimeLifted<<std::endl;
//	paramsFile<<"dense time gap lifted "<<denseTimeLifted<<std::endl;

//	if(parameters.count("NEGATIVE_THRESHOLD_LIFTED")>0){
//		negativeThresholdLifted=std::stod(parameters["NEGATIVE_THRESHOLD_LIFTED"]);
//	}
//	else{
//		negativeThresholdLifted=-1;
//	}
//	std::cout<<"negative threshold lifted "<<negativeThresholdLifted<<std::endl;
//	paramsFile<<"negative threshold lifted "<<negativeThresholdLifted<<std::endl;

//	if(parameters.count("POSITIVE_THRESHOLD_LIFTED")>0){
//		positiveThresholdLifted=std::stod(parameters["POSITIVE_THRESHOLD_LIFTED"]);
//	}
//	else{
//		positiveThresholdLifted=1;
//	}
//	std::cout<<"positive threshold lifted "<<positiveThresholdLifted<<std::endl;
//	paramsFile<<"positive threshold lifted "<<positiveThresholdLifted<<std::endl;

//	if(parameters.count("LONGER_LIFTED_INTERVAL")>0){
//		longerIntervalLifted=std::stoul(parameters["LONGER_LIFTED_INTERVAL"]);
//	}
//	else{
//		longerIntervalLifted=4;
//	}
//	std::cout<<"longer interval lifted "<<longerIntervalLifted<<std::endl;
//	paramsFile<<"longer interval lifted "<<longerIntervalLifted<<std::endl;


//	if(parameters.count("MAX_TIMEGAP_COMPLETE")>0){
//		maxTimeGapComplete=std::stoul(parameters["MAX_TIMEGAP_COMPLETE"]);
//	}
//	else{
//		maxTimeGapComplete=maxTimeFrame;
//	}
//	std::cout<<"max time gap complete "<<maxTimeGapComplete<<std::endl;
//	paramsFile<<"max time gap complete "<<maxTimeGapComplete<<std::endl;


//	paramsFile.close();




//}









}//End namespace
}//End namespace

#endif /* INCLUDE_DISJOINT_PATHS_INPUTCONFIG_HXX_ */
