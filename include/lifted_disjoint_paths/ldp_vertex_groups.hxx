#ifndef LDP_VERTEX_GROUPS_HXX
#define LDP_VERTEX_GROUPS_HXX

#include <stdexcept>
//#include "disjoint-paths/disjointPathsMethods.hxx"
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <andres/graph/digraph.hxx>
#include <iterator>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <utility>
#include <unordered_map>
#include "ldp_file_processing_methods.hxx"


namespace LPMP {


template<class T=size_t>
class VertexGroups{
public:
    VertexGroups(std::unordered_map<size_t,std::vector<size_t>> groups_,std::vector<size_t> vToGroup_):
        vToGroup(vToGroup_)
    {
        maxVertex=vToGroup_.size()-3;
        maxTime=*(vToGroup.rbegin())-1;
        groups=std::vector<std::vector<size_t>>(maxTime+2);
        for(auto pair:groups_){
            groups[pair.first]=pair.second;
        }

    }
    VertexGroups(){
        maxVertex=0;
        maxTime=0;
    }




    template<class PAR> VertexGroups(PAR& parameters){
        initFromFile(parameters.getTimeFileName(),parameters);

    }

    template<class PAR>
    void initFromFile(const std::string& fileName,const PAR& parameters);

    void initFromVector(const std::vector<size_t>& verticesInFrames);




    size_t getGroupIndex(size_t v) const{
        return vToGroup[v];
    }

    const std::vector<size_t>& getGroupVertices(size_t index)const{
        return groups.at(index);
    }

    //ID of maximal valid vertex (i.e. without s and t)
    size_t getMaxVertex() const {
        return maxVertex;
    }

    //Time of the last video frame
    size_t getMaxTime() const {
        return maxTime;
    }

    std::vector<std::vector<size_t>> extractInnerPaths(const std::vector<std::vector<size_t> > &paths, const size_t minT,const size_t maxT) const;


private:
    //std::vector<std::vector<size_t>> groups;
    std::vector<std::vector<size_t>> groups;
    std::vector<size_t> vToGroup;
    size_t maxVertex;
    size_t maxTime;


};


template<class T>
inline void VertexGroups<T>::initFromVector(const std::vector<size_t>& verticesInFrames){
    maxTime=verticesInFrames.size();
    groups=std::vector<std::vector<size_t>>(maxTime+2);
    size_t inFrameCounter=0;
    size_t vertexCounter=0;
    size_t frameCounter=1;
    std::vector<size_t> verticesInGroup;
    vToGroup=std::vector<size_t>();
    while(frameCounter<=maxTime){
        while(frameCounter<=maxTime&&inFrameCounter==verticesInFrames.at(frameCounter-1)){
            groups.at(frameCounter)=verticesInGroup;
            inFrameCounter=0;
            frameCounter++;
            verticesInGroup=std::vector<size_t>();
        }
        if(frameCounter<=maxTime){
            verticesInGroup.push_back(vertexCounter);
            vToGroup.push_back(frameCounter);
            inFrameCounter++;
            vertexCounter++;
        }
    }

    maxVertex=vertexCounter-1;
    size_t s=maxVertex+1;
    size_t t=maxVertex+2;

    verticesInGroup=std::vector<size_t>();
    verticesInGroup.push_back(s);
    vToGroup.push_back(0);
    groups.at(0)=verticesInGroup;

    verticesInGroup=std::vector<size_t>();
    verticesInGroup.push_back(t);
    vToGroup.push_back(frameCounter);
    groups.at(frameCounter)=verticesInGroup;

    for (int i = 0; i < groups.size(); ++i) {
        for (int j = 0; j < groups[i].size(); ++j) {
            assert(vToGroup.at(groups[i][j])==i);
        }
    }
    //    std::cout<<"max vertex "<<maxVertex<<std::endl;
    //    std::cout<<"max time "<<maxTime<<std::endl;



}


template<class T>
template<class PAR>
inline void VertexGroups<T>::initFromFile(const std::string& fileName, const PAR &parameters){
    size_t lineCounter=0;
    std::vector<size_t> currentGroup;
    std::vector<std::string> strings;
    std::string line;
    char delim=',';


    currentGroup=std::vector<size_t>();
    groups.push_back(currentGroup);
    currentGroup=std::vector<size_t>();


    size_t maxTimeToRead=parameters.getMaxTimeFrame();



    std::ifstream timeData;
    try{
        timeData.open(fileName);
        if(!timeData){
            throw std::system_error(errno, std::system_category(), "failed to open file with vertices in time layers "+fileName);
        }

        unsigned int previousTime=1;
        unsigned int time;

        while (std::getline(timeData, line) && !line.empty()) {


            lineCounter++;

            strings = split(line,delim);

            if (strings.size() < 2) {
                throw std::runtime_error(
                            std::string("Vertex and time frame expected"));
            }

            unsigned int v = std::stoul(strings[0]);
            time = std::stoul(strings[1]);
            if(time>maxTimeToRead){

                break;
            }
            if(vToGroup.size()!=v){
                throw std::runtime_error(
                            std::string("Wrong vertex numbering in time file"));
            }
            else{

                vToGroup.push_back(time);

                if(time==previousTime){
                    currentGroup.push_back(v);
                }

                else{

                    groups.push_back(currentGroup);
                    currentGroup=std::vector<size_t>();
                    while(groups.size()<time){
                        groups.push_back(currentGroup);
                        currentGroup=std::vector<size_t>();
                    }
                    currentGroup.push_back(v);
                }
                previousTime=time;


            }
        }
        groups.push_back(currentGroup);


        maxTime=*(vToGroup.rbegin());


        vToGroup.push_back(0);

        groups[0].push_back(vToGroup.size()-1);

        //time frame of t is maxTime
        vToGroup.push_back(maxTime+1);
        currentGroup=std::vector<size_t>();
        currentGroup.push_back(vToGroup.size()-1);
        groups.push_back(currentGroup);


        maxVertex=vToGroup.size()-3;

    }
    catch (std::system_error& er) {
        std::clog << er.what() << " (" << er.code() << ")" << std::endl;

    }

}

template<class T>
inline std::vector<std::vector<size_t>> VertexGroups<T>::extractInnerPaths(const std::vector<std::vector<size_t>>& paths, const size_t minT, const size_t maxT) const{
    //maxT inclusive
    std::vector<std::vector<size_t>> outputPaths;
    for (int i = 0; i < paths.size(); ++i) {
        const std::vector<size_t>& path=paths[i];
        std::vector<size_t> outputPath;
        for(size_t vertex: path){
            size_t time=getGroupIndex(vertex);
            if(time<=maxT){
                if(time>=minT){
                    outputPath.push_back(vertex);
                }
            }
            else break;
        }
        if(outputPath.size()>0){
            outputPaths.push_back(outputPath);
        }
    }
    return outputPaths;
}






}

#endif // LDP_VERTEX_GROUPS_HXX
