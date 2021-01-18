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
#include "ldp_parameter_parser.hxx"
#include <config.hxx>

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
class LdpParameters {
public:

//	~ConfigDisjoint(){
//        controlOutput<<"config disjoint destructor"<<std::endl;
//        writeControlOutput();
//        //infoFile()<<"config disjoint destructor"<<std::endl;
//        //(*pInfoFile).close();
//        //delete pInfoFile;
//	}


    LdpParameters(const std::string& inputFileName);
    LdpParameters(std::map<std::string,std::string>& parameters);


	const std::string& getGraphFileName() const {
		return graphFileName;
	}

	size_t getMaxTimeBase() const {
		return maxTimeBase;
	}

	size_t getMaxTimeFrame() const {
		return maxTimeFrame;
	}

    size_t getMinTimeFrame() const {
        return minTimeFrame;
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

    bool isUseAdaptiveThreshold() const{
        return useAdaptiveThresholds;
    }


    const double& getMustCutPenalty()const{
        return mustCutPenalty;
    }


//	void setMaxTimeLifted(size_t maxTimeLifted) {
//		this->maxTimeLifted = maxTimeLifted;
//	}


	size_t getMaxTimeGapComplete() const {
        return maxTimeGapComplete;
	}


//	std::ofstream& infoFile() const {
//		std::ofstream & infoF=*pInfoFile;
//		return infoF;
//	}

    std::stringstream& getControlOutput(){
          return controlOutput;
      }


    void writeControlOutput(){
          if(diagnostics()) std::cout<<controlOutput.str();

//            if(controlOutputFiles){
//                infoFile()<<controlOutput.str();
//                infoFile().flush();
//            }
            controlOutput=std::stringstream();
        }


	//ConfigDisjoint<>& operator=(const ConfigDisjoint<>&);


    size_t getTightenMaxEdgeUsage() const{
        return tighteningMaxEdgeUsage;
    }

    double getTightenMinImprovement() const{
        return tighteningMinImprovement;
    }


    bool isAllBaseZero()const {
        return allBaseToZero;
    }

    bool isBaseCoverdWithLifted()const{
        return coverBaseWithLifted;
    }

    bool isMustCutMissing()const{
        return missingAsMustCut;
    }


private:
    LdpParameters<T>(const LdpParameters<T>&);
    LdpParameters();

    void init(std::map<std::string,std::string>& parameters);

    //std::pair<std::string,std::string> parseLine(std::string line,char delim);
	//std::ofstream * fileForGeneralOutputs;
   // mutable std::ofstream* pInfoFile;
    std::string inputFileName;  //file with 3 paths: graph, time information, parameters

    std::string graphFileName;  //file with input graph structure
	std::string timeFileName;  //file with input time frames information

    std::string inputParamsFileName; //file where input parameters are read
	bool debugOutputFiles;  //use debugging output files

	double inputCost;   //cost of input edges (starting a track)
	double outputCost;  //cost of output edges (terminating a track)

	bool sparsify;     //use graph sparsification?
	bool restrictFrames;  //is maximal number of frames given
	size_t maxTimeFrame;   //maximal time frame to be used
    size_t minTimeFrame;
	//size_t smallIntervals;

    bool useAdaptiveThresholds;

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


     double tighteningMinImprovement;
     size_t tighteningMaxEdgeUsage;


     bool allBaseToZero;
     bool coverBaseWithLifted;

     bool missingAsMustCut;
     double mustCutPenalty;

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
inline LdpParameters<T>::LdpParameters(const std::string &inputFileName){
    ParametersParser parser;
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
inline LdpParameters<T>::LdpParameters(std::map<std::string,std::string>& parameters){
    init(parameters);
}

template<class T>
inline void LdpParameters<T>::init(std::map<std::string,std::string>& parameters){


    if(parameters.count("INPUT_FRAMES")>0){
        timeFileName=parameters["INPUT_FRAMES"];
    }
    else{
        timeFileName="problemDesc_frames";
    }
    controlOutput<<"input frames "<<timeFileName<<std::endl;
    writeControlOutput();
    //paramsFile<<"input frames "<<timeFileName<<std::endl;


    if(parameters.count("INPUT_GRAPH")>0){
        graphFileName=parameters["INPUT_GRAPH"];
    }
    else{
        graphFileName="";
        //throw std::runtime_error("File with the input graph was not specified in the config file");
    }
    controlOutput<<"input graph "<<graphFileName<<std::endl;

    writeControlOutput();




    if(parameters.count("SPARSIFY")>0){
        sparsify=std::stoi(parameters["SPARSIFY"]);
    }
    else{
        sparsify=1;
    }
    controlOutput<<"sparsify "<<sparsify<<std::endl;
    writeControlOutput();

    if(parameters.count("ALL_BASE_TO_ZERO")>0){
        allBaseToZero=std::stoi(parameters["ALL_BASE_TO_ZERO"]);
    }
    else{
        allBaseToZero=1;
    }
    controlOutput<<"all base to zero "<<allBaseToZero<<std::endl;
    writeControlOutput();


    if(parameters.count("COVER_BASE_WITH_LIFTED")>0){
        coverBaseWithLifted=std::stoi(parameters["COVER_BASE_WITH_LIFTED"]);
    }
    else{
        coverBaseWithLifted=1;
    }
    controlOutput<<"cover base with lifted "<<coverBaseWithLifted<<std::endl;
    writeControlOutput();




    if(parameters.count("ALL_BASE_TO_ZERO")>0){
        allBaseToZero=std::stoi(parameters["ALL_BASE_TO_ZERO"]);
    }
    else{
        allBaseToZero=1;
    }
    controlOutput<<"all base to zero "<<allBaseToZero<<std::endl;
    writeControlOutput();

    if(parameters.count("MAX_TIMEGAP")>0){
        maxTimeFrame=std::stoul(parameters["MAX_TIMEGAP"]);
        restrictFrames=1;
    }
    else{
        restrictFrames=0;
        maxTimeFrame=std::numeric_limits<size_t>::max();
    }
    controlOutput<<"max time frame "<<maxTimeFrame<<std::endl;
    writeControlOutput();

    if(parameters.count("MIN_TIME_FRAME")>0){
        minTimeFrame=std::stoul(parameters["MIN_TIME_FRAME"]);

    }
    else{

        minTimeFrame=1;
    }
    controlOutput<<"min time frame "<<minTimeFrame<<std::endl;
    writeControlOutput();


    if(parameters.count("INPUT_COST")>0){
        inputCost=std::stod(parameters["INPUT_COST"]);
    }
    else{
        inputCost=0;
    }
    controlOutput<<"input cost "<<inputCost<<std::endl;
    writeControlOutput();

    if(parameters.count("OUTPUT_COST")>0){
        outputCost=std::stod(parameters["OUTPUT_COST"]);
    }
    else{
        outputCost=0;
    }
    controlOutput<<"output cost "<<outputCost<<std::endl;
    writeControlOutput();


    if(parameters.count("MAX_TIMEGAP_BASE")>0){
        maxTimeBase=std::stoul(parameters["MAX_TIMEGAP_BASE"]);
    }
    else{
        maxTimeBase=60;
    }
    controlOutput<<"max time gap base "<<maxTimeBase<<std::endl;
    writeControlOutput();

    if(parameters.count("KNN_GAP")>0){
        knnTimeGap=std::stoul(parameters["KNN_GAP"]);
    }
    else{
        knnTimeGap=3;
    }
    controlOutput<<"KNN time gap "<<knnTimeGap<<std::endl;
    writeControlOutput();


    if(parameters.count("KNN_K")>0){
        knnK=std::stoul(parameters["KNN_K"]);
    }
    else{
        knnK=3;
    }
    controlOutput<<"KNN k "<<knnK<<std::endl;
    writeControlOutput();


    if(parameters.count("BASE_THRESHOLD")>0){
        baseUpperThreshold=std::stod(parameters["BASE_THRESHOLD"]);
    }
    else{
        baseUpperThreshold=0;
    }
    controlOutput<<"base upper threshold "<<baseUpperThreshold<<std::endl;
    writeControlOutput();


    if(parameters.count("MAX_TIMEGAP_LIFTED")>0){
        maxTimeLifted=std::stoul(parameters["MAX_TIMEGAP_LIFTED"]);
    }
    else{
        maxTimeLifted=60;
    }
    controlOutput<<"max time gap lifted "<<maxTimeLifted<<std::endl;
    writeControlOutput();

    if(parameters.count("DENSE_TIMEGAP_LIFTED")>0){
        denseTimeLifted=std::stoul(parameters["DENSE_TIMEGAP_LIFTED"]);
    }
    else{
        denseTimeLifted=20;
    }
    controlOutput<<"dense time gap lifted "<<denseTimeLifted<<std::endl;
    writeControlOutput();

    if(parameters.count("NEGATIVE_THRESHOLD_LIFTED")>0){
        negativeThresholdLifted=std::stod(parameters["NEGATIVE_THRESHOLD_LIFTED"]);
    }
    else{
        negativeThresholdLifted=-1;
    }
    controlOutput<<"negative threshold lifted "<<negativeThresholdLifted<<std::endl;
    writeControlOutput();

    if(parameters.count("POSITIVE_THRESHOLD_LIFTED")>0){
        positiveThresholdLifted=std::stod(parameters["POSITIVE_THRESHOLD_LIFTED"]);
    }
    else{
        positiveThresholdLifted=1;
    }
    controlOutput<<"positive threshold lifted "<<positiveThresholdLifted<<std::endl;
    writeControlOutput();


    if(parameters.count("LONGER_LIFTED_INTERVAL")>0){
        longerIntervalLifted=std::stoul(parameters["LONGER_LIFTED_INTERVAL"]);
    }
    else{
        longerIntervalLifted=4;
    }
    controlOutput<<"longer interval lifted "<<longerIntervalLifted<<std::endl;
    writeControlOutput();



    maxTimeGapComplete=0;
    if(parameters.count("MAX_TIMEGAP_COMPLETE")>0){
        maxTimeGapComplete=std::stoul(parameters["MAX_TIMEGAP_COMPLETE"]);
    }
    maxTimeGapComplete=std::max(maxTimeGapComplete,maxTimeLifted);
    maxTimeGapComplete=std::max(maxTimeGapComplete,maxTimeBase);
//    else{
//        maxTimeGapComplete=maxTimeFrame;
//    }
    controlOutput<<"max time gap complete "<<maxTimeGapComplete<<std::endl;

    if(parameters.count("USE_ADAPTIVE_THRESHOLDS")){
        useAdaptiveThresholds=std::stoi(parameters["USE_ADAPTIVE_THRESHOLDS"]);
    }
    else {
        useAdaptiveThresholds=false;
    }
     controlOutput<<"adaptive thresholds "<<useAdaptiveThresholds<<std::endl;


     if(parameters.count("TIGHT_MIN_IMPROVEMENT")>0){
         tighteningMinImprovement=std::stod(parameters["TIGHT_MIN_IMPROVEMENT"]);
     }
     else{
         tighteningMinImprovement=0.00001;
     }
     assert(tighteningMinImprovement>=0);
     controlOutput<<"minimal improvement for tightening "<<tighteningMinImprovement<<std::endl;

     if(parameters.count("TIGHT_MAX_EDGE_USAGE")>0){
         tighteningMaxEdgeUsage=std::stoi(parameters["TIGHT_MAX_EDGE_USAGE"]);
     }
     else{
         tighteningMaxEdgeUsage=4;
     }
     controlOutput<<"maximal edge usage for tightening "<<tighteningMaxEdgeUsage<<std::endl;

     if(parameters.count("MISSING_AS_MUST_CUT")>0){
         missingAsMustCut=std::stoi(parameters["MISSING_AS_MUST_CUT"]);
     }
     else{
         missingAsMustCut=0;
     }
     controlOutput<<"missing edges as must cut "<<missingAsMustCut<<std::endl;


     if(parameters.count("MUST_CUT_PENALTY")>0){
         mustCutPenalty=std::stod(parameters["MUST_CUT_PENALTY"]);
     }
     else{
         mustCutPenalty=100.0;
     }
     controlOutput<<"must cut penalty "<<mustCutPenalty<<std::endl;


     writeControlOutput();


}







}//End namespace
}//End namespace

#endif /* INCLUDE_DISJOINT_PATHS_INPUTCONFIG_HXX_ */
