#ifndef LDP_INTERVAL_CONNECTION_HXX
#define LDP_INTERVAL_CONNECTION_HXX

#include "ldp_directed_graph.hxx"
#include "ldp_vertex_groups.hxx"
#include "ldp_paths_extractor.hxx"
#include "ldp_functions.hxx"
#include <chrono>


namespace LPMP {

class LdpLabelsAssignment{
public:
    LdpLabelsAssignment(){
        maxUsedLabel=0;
    }

    const std::vector<std::array<size_t,2>> getVerticesToLabels()const{
        return verticesToLabels;
    }
    void init(const std::vector<std::vector<size_t>>& initialPaths){
        assert(maxUsedLabel==0);
        assert(lastPathsLabels.size()==0);
        assert(verticesToLabels.size()==0);
        lastPathsLabels=std::vector<size_t>(initialPaths.size());
        for (size_t i = 0; i < initialPaths.size(); ++i) {
            for (size_t j = 0; j < initialPaths[i].size(); ++j) {
                size_t vertex=initialPaths[i][j];
                size_t label=i+1;
                verticesToLabels.push_back({vertex,label});

            }
            lastPathsLabels[i]=i+1;
        }
        maxUsedLabel=initialPaths.size();
    }

    void update(const std::vector<std::vector<size_t>>& newPaths,const std::vector<int>& oldToNewPaths){
        assert(oldToNewPaths.size()==lastPathsLabels.size());
        std::vector<size_t> newPathsLabels(newPaths.size());

        std::vector<bool> processedPaths(newPaths.size(),false);

        for (size_t i = 0; i < oldToNewPaths.size(); ++i) {
            if(oldToNewPaths[i]!=-1){
                int index=oldToNewPaths[i];
                assert(index>=0&&index<newPaths.size());
                size_t label=lastPathsLabels[i];
                newPathsLabels[index]=label;
                processedPaths[index]=true;
                for (size_t j = 0; j < newPaths[index].size(); ++j) {
                    size_t vertex=newPaths[index][j];
                    verticesToLabels.push_back({vertex,label});

                }

            }
        }

        for (size_t i = 0; i < newPaths.size(); ++i) {
            if(!processedPaths[i]){
                maxUsedLabel++;
                size_t label=maxUsedLabel;
                for (size_t j = 0; j < newPaths[i].size(); ++j) {
                    size_t vertex=newPaths[i][j];
                    verticesToLabels.push_back({vertex,label});
                }
            }
        }
        lastPathsLabels=newPathsLabels;

    }

private:
    std::vector<std::array<size_t,2>> verticesToLabels;
    std::vector<size_t> lastPathsLabels;
    size_t maxUsedLabel;
};



//class LdpIntervalConnectionResult{
//public:
//    LdpIntervalConnectionResult(const std::vector<std::vector<size_t>>& _firstIntervalPaths,
//    const std::vector<std::vector<size_t>>& _secondIntervalPaths):
//        firstIntervalPaths(_firstIntervalPaths),
//        secondIntervalPaths(_secondIntervalPaths)
//    {
//    }

//    const std::vector<std::vector<size_t>>& getFirstIntervalPaths()const{
//        return firstIntervalPaths;
//    }

//    const std::vector<std::vector<size_t>>& getSecondIntervalPaths()const{
//        return secondIntervalPaths;
//    }

//    void setMiddleIntervalPaths(const std::vector<std::vector<size_t>>& paths){
//        middleIntervalPaths=paths;
//    }

//    const std::vector<std::vector<size_t>>& getMiddleIntervalPaths()const{
//        return middleIntervalPaths;
//    }

//    void setFirstToMiddle(const std::vector<int>& ftm){
//        firstToMiddle=ftm;
//    }

//    const std::vector<int>& getFirstToMiddle()const{
//        return firstToMiddle;
//    }

//    void setMiddleToSecond(const std::vector<int>& mts){
//        middleToSecond=mts;
//    }

//    const std::vector<int>& getMiddleToSecond()const{
//        return middleToSecond;
//    }


//    void setStartInMiddle(const std::vector<size_t>& sim){
//        startInMiddle=sim;
//    }


//    const std::vector<size_t>& getStartInMiddle()const{
//        return startInMiddle;
//    }

//    void setStartInSecond(const std::vector<size_t>& sis){
//        startInSecond=sis;
//    }

//    const std::vector<size_t>& getStartInSecond()const{
//        return startInSecond;
//    }



//private:
//    const std::vector<std::vector<size_t>>& firstIntervalPaths;
//    const std::vector<std::vector<size_t>>& secondIntervalPaths;
//    std::vector<std::vector<size_t>> middleIntervalPaths;
//    std::vector<int> firstToMiddle;
//    std::vector<int> middleToSecond;
//    std::vector<size_t> startInMiddle;
//    std::vector<size_t> startInSecond;

//};


class LdpIntervalConnection{
public:
     LdpIntervalConnection(  const LdpPathsExtractor& _pathExtractor1,
                            const LdpPathsExtractor& _pathExtractor2);

     void initEdgesFromVectors(const std::vector<std::array<size_t,2> > &edges, const std::vector<double> &costs);

     void initFromFile(const std::string& fileName);

     void initScoreOfVertices(const std::vector<size_t> &listOfVertices, const std::vector<double> &costs);

     const LdpDirectedGraph& getCompleteGraph()const{
         return completeGraph;
     }
     const VertexGroups<>& getVertexGroups()const{
         return vg;
     }

     const std::vector<double> getScoreOfVertices()const{
         return verticesScore;
     }

     const std::chrono::steady_clock::time_point& getContructorBegin()const {
         return constructorBegin;
     }

     std::vector<std::vector<size_t>> decodePaths(const std::vector<std::vector<size_t> > &paths)const;

     void createResultsStructures(const std::vector<std::vector<size_t>> &paths);


     const std::vector<std::vector<size_t>>& getFirstIntervalPaths()const{
         return pathExtractor1.getExtractedPaths();
     }

     const std::vector<std::vector<size_t>>& getSecondIntervalPaths()const{
         return pathExtractor2.getExtractedPaths();
     }


     const std::vector<std::vector<size_t>>& getMiddleIntervalPaths()const{
         return middleIntervalPaths;
     }


     const std::vector<int>& getFirstToMiddle()const{
         return firstToMiddle;
     }

     const std::vector<int>& getMiddleToSecond()const{
         return middleToSecond;
     }



private:

    size_t localToGlobalFree(const size_t& index)const {
        assert(index<n1+n2);
        assert(index>=n1);
        return index-n1+firstFreeVertex;

    }

    size_t globalToLocalFree(const size_t& id)const{
        assert(id>=firstFreeVertex&&id<=lastFreeVertex);
        return id-firstFreeVertex+n1;
    }

    size_t firstPathsVertexToLocalIndex(const size_t& vertexGlobalID)const{
        return pathExtractor1.pathToVertex(vertexGlobalID);
    }

    size_t secondPathsVertexToLocalIndex(const size_t& vertexGlobalID)const{
        assert(pathExtractor2.pathToVertex(vertexGlobalID)+n1+n2<numberOfLocalVertices);
        return pathExtractor2.pathToVertex(vertexGlobalID)+n1+n2;
    }

    size_t secondPathsLocalIndexToPathIndex(const size_t& vertexIndex)const{
        assert(vertexIndex>=n1+n2&&vertexIndex<numberOfLocalVertices);
        return vertexIndex-n1-n2;
    }


   char vertexToGraphPartGlobal(const size_t& vertexGlobalID)const;
   char vertexToGraphPartLocal(const size_t& localVertexIndex) const;



    const LdpPathsExtractor& pathExtractor1;
    const LdpPathsExtractor& pathExtractor2;
    std::vector<double> verticesScore;
    size_t firstFreeVertex;
    size_t lastFreeVertex;
    //const std::vector<std::vector<size_t>>& firstPaths;
    //const std::vector<std::vector<size_t>>& secondPaths;
    //std::vector<size_t> vToPathFirst;  //Shifted by firstMinVertex.  What is the number of the fixed path which the vertex belongs to
   // std::vector<size_t> vToPathSecond; //Shifted by secondMinVertex.
    size_t n1;
    size_t n2;
    size_t n3;
    size_t numberOfLocalVertices;
    LdpDirectedGraph completeGraph;
    VertexGroups<> vg;
    std::chrono::steady_clock::time_point constructorBegin;
    size_t minVertex;



   // const std::vector<std::vector<size_t>>& firstIntervalPaths;
    //const std::vector<std::vector<size_t>>& secondIntervalPaths;
    std::vector<std::vector<size_t>> middleIntervalPaths;
    std::vector<int> firstToMiddle;
    std::vector<int> middleToSecond;




};





//    size_t minTime;
//    size_t maxTime;
//    size_t minVertex;



//    minTime=1;
//    assert(vgInterval1.getMaxTime()>cutoffLength);
//    maxTime=vgInterval1.getMaxTime()-cutoffLength;
//    minVertex=firstMinVertex;
//    firstPaths=vgInterval1.extractInnerPaths(firstPaths,minTime,maxTime,firstMinVertex);


//    minTime=cutoffLength+1;
//    maxTime=vgInterval2.getMaxTime();
//    assert(minTime<maxTime);
//    minVertex=secondMinVertex;
//    secondPaths=vgInterval2.extractInnerPaths(secondInputPahts,minTime,maxTime,secondMinVertex);


}

#endif // LDP_INTERVAL_CONNECTION_HXX
