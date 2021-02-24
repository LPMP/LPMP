#ifndef LDP_PRIMAL_HEURISTICS_HXX
#define LDP_PRIMAL_HEURISTICS_HXX

#include <cstdlib>
#include <vector>
#include <lifted_disjoint_paths/ldp_instance.hxx>
#include "graph_matching/graph_matching_input.h"
#include "graph_matching/min_cost_flow_factor_ssp.hxx"
#include "graph_matching/matching_problem_input.h"

namespace LPMP{

struct LdpHeuristicConnections{
public:
    size_t pathIndex1;
    size_t pathIndex2;
  //  size_t alternativeFirstEnd;  //Not needed, the cuts will be done before the merging starts
  //  size_t alternativeSecondBegin;
    double connectionValue;

};

template <class SNC_FACTOR>
class LdpPrimalHeuristics{
public:
    LdpPrimalHeuristics(const std::vector<size_t>& _vertexLabels,
                        const std::vector<size_t>& _startingVertices,
                        //const std::vector<std::vector<size_t>>& _paths,
                        const std::vector<size_t>& _neighboringVertices,
                        const lifted_disjoint_paths::LdpInstance* _pInstance,
                        std::vector<std::array<SNC_FACTOR*,2>>* _p_single_node_cut_factors,
                        bool useRepamCost=false);

    ~LdpPrimalHeuristics(){
        if(useReparametrizedCosts) delete pRepamLiftedGraph;
    }



    const std::vector<size_t>& getVertexLabels()const{
        return vertexLabels;
    }

    const std::vector<std::vector<size_t>>& getPaths()const{
        return adjustedPaths;
    }

    const std::vector<size_t>& getStartingVertices()const{
        return startingVertices;
    }

    const std::vector<size_t>& getNeighboringVertices()const{
        return neighboringVertices;
    }

    double getPrimalValue()const{
        return currentPrimalValue;
    }





    void evaluateAll();

void improveAllPaths();


private:

    std::array<double,2> mutualCostToSubtract(size_t vertex,size_t otherPathVertex,size_t targetPathIndex, bool useForward);
    void finalizeResults(bool changeSNC);

    double currentPrimal(bool checkImprovement=true);
    void cutOnePathWithMutualCostUpdate(size_t cutVertex);


    void initReverseNeighbors();

    std::array<size_t,2> findBestCut(size_t pathIndex1,size_t pathIndex2);
    void cutToEnableConnections();
   // void cutToEnableConnectionsNew();
    void initLastVertices();
    void mergePaths();

    void connectToOnePath(size_t indexOfPath, const std::vector<std::map<size_t,double>>& costsBetweenPaths, size_t bestNeighborIndex);

    void computeCummulativeCostsGlobal();
    void initCutValues();
    void updateCutValues(const std::vector<size_t>& oldLabels, const std::vector<size_t>& oldNeighbors);
    std::vector<size_t> findCutCandidates(const std::vector<double>& cutValues);
    void cutPaths(const std::vector<size_t>& cutCandidates,const std::vector<size_t> indicesOfPathsToCut);

    void initCummulativeCosts();

    const LdpDirectedGraph& getLiftedGraph(){
        if(!useReparametrizedCosts){
            return pInstance->getMyGraphLifted();
        }
        else{
            return (*pRepamLiftedGraph);
        }
    }




    double getInputCost(size_t nodeID){
        assert(nodeID<pInstance->getNumberOfVertices()-2);
        const size_t& s=pInstance->getSourceNode();
        assert(baseEdges.size()>s);
        assert(baseEdges[s].count(nodeID)>0);
        return baseEdges[s][nodeID];

    }

    double getOutputCost(size_t nodeID){
        assert(nodeID<baseEdges.size());
        const size_t& t=pInstance->getTerminalNode();
        auto it=baseEdges[nodeID].rbegin();
        assert(it->first==t);
        return it->second;
    }



    size_t maxMoveForEndCut;

    std::vector<std::map<size_t,double>> baseEdges;
    std::vector<size_t> vertexLabels;
    std::vector<size_t> startingVertices;
    std::vector<std::vector<size_t>> adjustedPaths;
    std::vector<size_t> neighboringVertices;
    std::vector<size_t> reverseNeighbors;
    const lifted_disjoint_paths::LdpInstance* pInstance;
    size_t numberOfPaths;
    double currentPrimalValue;

    std::vector<size_t> lastVertices;

    std::vector<std::array<SNC_FACTOR*,2>>* pSNCFactors;

    size_t currentTime;


     std::vector<size_t> pointersToPaths;
    std::vector<std::vector<double>> cummulativeCosts;  //dimension number of paths x number of paths, really the full dimension!
    std::vector<std::vector<double>> negativeMutualCosts;
    std::vector<double> cutValues;

    std::vector<size_t> checkArray;
    std::map<size_t,std::set<size_t>> checkMap;

    bool useReparametrizedCosts;
    LdpDirectedGraph * pRepamLiftedGraph;
    std::vector<double> costsOfPaths;

    double mergeThreshold;
    //LdpDirectedGraph * pRepamBaseGraph;

};


template<class SNC_FACTOR>
LdpPrimalHeuristics<SNC_FACTOR>::LdpPrimalHeuristics(const std::vector<size_t>& _vertexLabels,
                    const std::vector<size_t>& _startingVertices,
                    //const std::vector<std::vector<size_t>>& _paths,
                    const std::vector<size_t>& _neighboringVertices,
                    const lifted_disjoint_paths::LdpInstance* _pInstance,
                    std::vector<std::array<SNC_FACTOR*,2>>* _p_single_node_cut_factors,
                    bool useRepamCost)
  {
    useReparametrizedCosts=useRepamCost;
    vertexLabels=_vertexLabels;
    startingVertices=_startingVertices;
   // paths=_paths;
    neighboringVertices=_neighboringVertices;
    pInstance=_pInstance;
    pSNCFactors=_p_single_node_cut_factors;

    numberOfPaths=startingVertices.size();
    const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();

    baseEdges=std::vector<std::map<size_t,double>>(baseGraph.getNumberOfVertices());

    currentPrimalValue=0;

    if(useReparametrizedCosts){
       pRepamLiftedGraph=new LdpDirectedGraph(pInstance->getMyGraphLifted()); //TODO new, copy structure from orig graphs, set costs
       //pRepamBaseGraph=new LdpDirectedGraph(pInstance->getMyGraph());  //TODO: probably can be removed and use baseEdges vector of maps
       for (size_t i = 0; i < pRepamLiftedGraph->getNumberOfVertices(); ++i) {
           size_t numberOfNeighbors=pRepamLiftedGraph->getNumberOfEdgesFromVertex(i);
           auto pOutFactor=(*pSNCFactors)[i][1]->get_factor();
           const std::vector<double>& liftedCosts=pOutFactor->getLiftedCosts();
           assert(liftedCosts.size()==numberOfNeighbors);
           for (size_t j = 0; j < numberOfNeighbors; ++j) {
               pRepamLiftedGraph->setForwardEdgeCost(i,j,liftedCosts[j]);
               assert(pRepamLiftedGraph->getForwardEdgeVertex(i,j)==pOutFactor->getLiftedIDs()[j]);
           }
       }
       for (size_t i = 0; i < pRepamLiftedGraph->getNumberOfVertices(); ++i) {
           size_t numberOfNeighbors=pRepamLiftedGraph->getNumberOfEdgesToVertex(i);
           auto pInFactor=(*pSNCFactors)[i][0]->get_factor();
           const std::vector<double>& liftedCosts=pInFactor->getLiftedCosts();
           assert(liftedCosts.size()==numberOfNeighbors);
           for (size_t j = 0; j < numberOfNeighbors; ++j) {
               pRepamLiftedGraph->updateBackwardEdgeCost(i,j,liftedCosts[j]);
               double controlCost=pInstance->getMyGraphLifted().getBackwardEdgeCost(i,j);
             //  assert(abs(controlCost-pRepamLiftedGraph->getBackwardEdgeCost(i,j))<eps);
               assert(pRepamLiftedGraph->getBackwardEdgeVertex(i,j)==pInFactor->getLiftedIDs()[j]);
           }
       }

       for (size_t i = 0; i < baseGraph.getNumberOfVertices()-2; ++i) {
           //size_t numberOfNeighbors=pRepamLiftedGraph->getNumberOfEdgesFromVertex(i);
           auto pOutFactor=(*pSNCFactors)[i][1]->get_factor();
           double nodeCost=pOutFactor->getNodeCost();
           const std::vector<double>& baseCosts=pOutFactor->getBaseCosts();
           const std::vector<size_t>& baseIDs=pOutFactor->getBaseIDs();
           size_t numberOfNeighbors=baseCosts.size();
           for (size_t j = 0; j < numberOfNeighbors; ++j) {
               baseEdges[i][baseIDs[j]]=baseCosts[j];
               baseEdges[i][baseIDs[j]]+=nodeCost;

           }
       }

       for (size_t i = 0; i < baseGraph.getNumberOfVertices()-2; ++i) {
           //size_t numberOfNeighbors=pRepamLiftedGraph->getNumberOfEdgesFromVertex(i);
           auto pInFactor=(*pSNCFactors)[i][0]->get_factor();
           double nodeCost=pInFactor->getNodeCost();
           const std::vector<double>& baseCosts=pInFactor->getBaseCosts();
           const std::vector<size_t>& baseIDs=pInFactor->getBaseIDs();
           size_t numberOfNeighbors=baseCosts.size();
           for (size_t j = 0; j < numberOfNeighbors; ++j) {
               baseEdges[baseIDs[j]][i]+=baseCosts[j];
               baseEdges[baseIDs[j]][i]+=nodeCost;
               double controlCost=pInstance->getMyGraph().getBackwardEdgeCost(i,j);
              // assert(abs(controlCost-baseEdges[baseIDs[j]][i])<eps);

           }
       }
       // pRepamBaseGraph;
    }
    else{
        pRepamLiftedGraph=nullptr;
        //pRepamBaseGraph=nullptr;

        for (int i = 0; i < baseGraph.getNumberOfVertices(); ++i) {
            for (auto iter=baseGraph.forwardNeighborsBegin(i);iter!=baseGraph.forwardNeighborsEnd(i);iter++) {
                baseEdges[i][iter->first]=iter->second;
            }

        }
    }


   // cutAllPaths();
    currentPrimal(false);
    if(diagnostics()) std::cout<<"init primal value "<<currentPrimalValue<<std::endl;


    currentTime=1;



   pointersToPaths=startingVertices;


   maxMoveForEndCut=10;

   mergeThreshold=pInstance->parameters.getMergeThreshold();





}

template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::initCummulativeCosts(){
    cummulativeCosts=std::vector<std::vector<double>>(numberOfPaths);
    negativeMutualCosts=std::vector<std::vector<double>>(numberOfPaths);
    for (int i = 0; i < numberOfPaths; ++i) {
        cummulativeCosts[i]=std::vector<double>(numberOfPaths);
        negativeMutualCosts[i]=std::vector<double>(numberOfPaths);
    }
}



template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::initReverseNeighbors(){
    reverseNeighbors=std::vector<size_t>(neighboringVertices.size(),pInstance->getSourceNode());
    size_t numberOfVertices=neighboringVertices.size();
    for (size_t i = 0; i < numberOfVertices; ++i) {
        if(vertexLabels[i]>0){
            size_t neighbor=neighboringVertices[i];
            if(neighbor<numberOfVertices){
                reverseNeighbors[neighbor]=i;
            }
        }
    }

    for (size_t i = 0; i < startingVertices.size(); ++i) {
        if(startingVertices[i]!=pInstance->getTerminalNode()){
            reverseNeighbors[startingVertices[i]]=pInstance->getSourceNode();
        }
    }
    //TODO for starting vertices neighbor is s

}


template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::computeCummulativeCostsGlobal(){
    initCummulativeCosts();

    const LdpDirectedGraph& liftedGraph=getLiftedGraph();


    //Addding new lifted edges to cummulative costs
    for (int i = 0; i < liftedGraph.getNumberOfVertices(); ++i) {
        size_t label=vertexLabels[i];
        if(label>0){
            auto iter=liftedGraph.forwardNeighborsBegin(i);
            for (;iter!=liftedGraph.forwardNeighborsEnd(i);iter++) {
                size_t vertex2=iter->first;

                size_t label2=vertexLabels[vertex2];

                if(label2>0){
                    assert(label2>0);
                    assert(label-1<cummulativeCosts.size());
                    assert(label2-1<cummulativeCosts[label-1].size());
                    cummulativeCosts[label-1][label2-1]+=iter->second;

                    if(iter->second<0){
                        negativeMutualCosts[label-1][label2-1]+=iter->second;
                    }
                }
            }
        }
    }

}



template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::connectToOnePath(size_t indexOfPath, const std::vector<std::map<size_t,double>>& costsBetweenPaths,size_t bestNeighborIndex){


    if(diagnostics()) std::cout<<"merging "<<indexOfPath<<", "<<bestNeighborIndex<<std::endl;

    const VertexGroups<>& vg=(pInstance->vertexGroups);

    size_t maxTimeToUse=vg.getGroupIndex(startingVertices[bestNeighborIndex]);
    std::vector<size_t> bestPathSoFar={indexOfPath};

    double negativeCost=negativeMutualCosts[indexOfPath][bestNeighborIndex];
    double positiveCost=cummulativeCosts[indexOfPath][bestNeighborIndex]-negativeMutualCosts[indexOfPath][bestNeighborIndex];


    bestPathSoFar.push_back(bestNeighborIndex);

    if(debug()){
        std::cout<<"merging path "<<indexOfPath<<" with "<<bestNeighborIndex<<" using "<<std::endl;
        for (int i = 0; i < bestPathSoFar.size(); ++i) {
            std::cout<<bestPathSoFar[i]<<",";
        }
        std::cout<<std::endl;
        std::cout<<"merge cost "<<costsBetweenPaths[indexOfPath].at(bestNeighborIndex)<<std::endl;
        std::cout<<"paths cummulative costs "<<cummulativeCosts[indexOfPath][indexOfPath]<<" "<<cummulativeCosts[bestNeighborIndex][bestNeighborIndex]<<std::endl;
    }

    assert(bestPathSoFar.front()==indexOfPath);

    size_t indexToConnect=1;
    size_t currentLastVertex=lastVertices[indexOfPath];
    size_t labelToAssign=indexOfPath+1;
    while(indexToConnect<bestPathSoFar.size()){
        assert(neighboringVertices[currentLastVertex]==pInstance->getTerminalNode());

        size_t pathToConnect=bestPathSoFar[indexToConnect];
        assert(pathToConnect<numberOfPaths);
        size_t vertexToConnect=startingVertices[pathToConnect];
        neighboringVertices[currentLastVertex]=vertexToConnect;
        startingVertices[pathToConnect]=pInstance->getTerminalNode();
        while(vertexToConnect!=pInstance->getTerminalNode()){
            vertexLabels[vertexToConnect]=labelToAssign;
            vertexToConnect=neighboringVertices[vertexToConnect];
        }
        currentLastVertex=lastVertices[pathToConnect];
        lastVertices[pathToConnect]=pInstance->getTerminalNode();
        lastVertices[indexOfPath]=currentLastVertex;
        indexToConnect++;

    }


    for (size_t i = 0; i < numberOfPaths; ++i) {
        for (size_t j = 1; j < bestPathSoFar.size(); ++j) {
            size_t oldPathIndex=bestPathSoFar[j];


            cummulativeCosts[i][indexOfPath]+=cummulativeCosts[i][oldPathIndex];
            negativeMutualCosts[i][indexOfPath]+=negativeMutualCosts[i][oldPathIndex];

            cummulativeCosts[indexOfPath][i]+=cummulativeCosts[oldPathIndex][i];
            negativeMutualCosts[indexOfPath][i]+=negativeMutualCosts[oldPathIndex][i];



        }

    }


    for (size_t i = 0; i < numberOfPaths; ++i) {
        for (size_t j = 1; j < bestPathSoFar.size(); ++j) {
            size_t oldPathIndex=bestPathSoFar[j];


            cummulativeCosts[i][oldPathIndex]=0;
            negativeMutualCosts[i][oldPathIndex]=0;

            cummulativeCosts[oldPathIndex][i]=0;
            negativeMutualCosts[oldPathIndex][i]=0;

        }

    }
}




template<class SNC_FACTOR>
std::array<double, 2> LdpPrimalHeuristics<SNC_FACTOR>::mutualCostToSubtract(size_t vertex,size_t otherPathVertex,size_t targetPathIndex, bool useForward){
    const LdpDirectedGraph& liftedGraph=getLiftedGraph();
    const LdpDirectedGraph::edge * beginPointer=nullptr;
    const LdpDirectedGraph::edge * endPointer=nullptr;

    size_t timeOfOtherVertex=pInstance->vertexGroups.getGroupIndex(otherPathVertex);
    if(useForward){
        beginPointer=liftedGraph.forwardNeighborsBegin(vertex);
        endPointer=liftedGraph.forwardNeighborsEnd(vertex);
      //  std::cout<<"forward"<<std::endl;
    }
    else{
        beginPointer=liftedGraph.backwardNeighborsBegin(vertex);
        endPointer=liftedGraph.backwardNeighborsEnd(vertex);
        // std::cout<<"backward"<<std::endl;
    }

    double valueToSubtract=0;
    double valueNegativeToSubtract=0;

    for (auto iter=beginPointer;iter!=endPointer;iter++) {
        size_t vertex2=iter->first;
        if(vertexLabels[vertex2]==targetPathIndex+1){
            size_t vertex2Time=pInstance->vertexGroups.getGroupIndex(vertex2);
            if(useForward&&vertex2Time>=timeOfOtherVertex){
              //  std::cout<<vertex<<","<<vertex2<<":"<<iter->second<<std::endl;
                valueToSubtract+=iter->second;
                if(iter->second<0){
                    valueNegativeToSubtract+=iter->second;
                }
            }
            else if(!useForward&&vertex2Time<=timeOfOtherVertex){
               // std::cout<<vertex<<","<<vertex2<<":"<<iter->second<<std::endl;
                valueToSubtract+=iter->second;
                if(iter->second<0){
                    valueNegativeToSubtract+=iter->second;
                }
            }
        }
    }

    return {valueToSubtract,valueNegativeToSubtract};
}

template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::cutOnePathWithMutualCostUpdate(size_t cutVertex){
    size_t oldPathIndex=vertexLabels[cutVertex];
    assert(oldPathIndex>0);
    oldPathIndex--;
    assert(cutVertex!=pInstance->getTerminalNode());

    size_t nextVertex=neighboringVertices[cutVertex];
    assert(nextVertex!=pInstance->getTerminalNode());

    assert(vertexLabels[nextVertex]==oldPathIndex+1);
    neighboringVertices[cutVertex]=pInstance->getTerminalNode();
    startingVertices.push_back(nextVertex);
    size_t newLabel=startingVertices.size();
    size_t newPathIndex=newLabel-1;

    while(nextVertex!=pInstance->getTerminalNode()){

        vertexLabels[nextVertex]=newLabel;
        nextVertex=neighboringVertices[nextVertex];
    }

    lastVertices.push_back(lastVertices[oldPathIndex]);
    lastVertices[oldPathIndex]=cutVertex;


    //TODO also update lastVertices
    const LdpDirectedGraph& liftedGraph=getLiftedGraph(); //maybe can be in a separate method

    for (size_t i = 0; i < numberOfPaths; ++i) {
        cummulativeCosts[i].push_back(0.0);
        negativeMutualCosts[i].push_back(0.0);
    }
    std::vector<double> newLine(numberOfPaths+1);
    cummulativeCosts.push_back(newLine);
    negativeMutualCosts.push_back(newLine);

    assert(cummulativeCosts.size()==newPathIndex+1);


    nextVertex=startingVertices.back();
    while(nextVertex!=pInstance->getTerminalNode()){
        for (auto iter=liftedGraph.forwardNeighborsBegin(nextVertex);iter!=liftedGraph.forwardNeighborsEnd(nextVertex);iter++) {
            size_t vertex2=iter->first;
            size_t label2=vertexLabels[vertex2];
            if(label2>0&&label2-1!=oldPathIndex){
                double cost=iter->second;

                cummulativeCosts[newPathIndex][label2-1]+=cost;
               if(label2-1!=newPathIndex) cummulativeCosts[oldPathIndex][label2-1]-=cost;
                if(cost<0){
                    negativeMutualCosts[newPathIndex][label2-1]+=cost;

                   if(label2-1!=newPathIndex) negativeMutualCosts[oldPathIndex][label2-1]-=cost;
                }
            }

        }
        for (auto iter=liftedGraph.backwardNeighborsBegin(nextVertex);iter!=liftedGraph.backwardNeighborsEnd(nextVertex);iter++) {
            size_t vertex2=iter->first;
            size_t label2=vertexLabels[vertex2];
            if(label2>0&&label2-1!=newPathIndex){
                double cost=iter->second;

                if(label2-1==oldPathIndex){
                    size_t v=vertex2;
                    while(v!=pInstance->getTerminalNode()){
                        cutValues[v]-=iter->second;
                        v=neighboringVertices[v];
                    }
                    v=startingVertices[newPathIndex];
                    while(v!=nextVertex){
                        cutValues[v]-=iter->second;
                        v=neighboringVertices[v];
                    }
                }


                cummulativeCosts[label2-1][newPathIndex]+=cost;
                cummulativeCosts[label2-1][oldPathIndex]-=cost;
                if(cost<0){
                    negativeMutualCosts[label2-1][newPathIndex]+=cost;

                    negativeMutualCosts[label2-1][oldPathIndex]-=cost;
                }
            }

        }
        nextVertex=neighboringVertices[nextVertex];
    }

    numberOfPaths=startingVertices.size();

}



template<class SNC_FACTOR>
std::array<size_t,2> LdpPrimalHeuristics<SNC_FACTOR>::findBestCut(size_t pathIndex1,size_t pathIndex2){
    size_t pointerFirst=lastVertices[pathIndex1];
    size_t pointerSecond=startingVertices[pathIndex2];
    assert(pointerFirst!=pInstance->getTerminalNode());
    assert(pointerSecond!=pInstance->getTerminalNode());

    double currentMutualCost=cummulativeCosts[pathIndex1][pathIndex2];
    double currentNegativeCost=negativeMutualCosts[pathIndex1][pathIndex2];



    double firstCutValue=0; //keep zeros even if including in/out costs
    double secondCutValue=0;
    size_t shiftsDone=0;

    std::array<size_t,2> baseEdgeVertices={pInstance->getTerminalNode(),pInstance->getTerminalNode()};

    size_t predFirst=reverseNeighbors[pointerFirst];
    size_t descSecond=neighboringVertices[pointerSecond];


    if(predFirst!=pInstance->getSourceNode()&&descSecond!=pInstance->getTerminalNode()){
        bool proceedSearch=true;


        std::array<double,2> lossOfFirstCost=mutualCostToSubtract(pointerFirst,pointerSecond,pathIndex2,true); //TODO take into account what has been removed on the other side
        std::array<double,2> lossOfSecondCost=mutualCostToSubtract(pointerSecond,pointerFirst,pathIndex1,false);


        double firstCutShifted=-cutValues[predFirst]-lossOfFirstCost[0]-baseEdges[predFirst].at(pointerFirst)+baseEdges[pInstance->getSourceNode()].at(pointerFirst);
        double secondCutShifted=-cutValues[pointerSecond]-lossOfSecondCost[0]-baseEdges[pointerSecond].at(descSecond)+baseEdges[pointerSecond].at(pInstance->getTerminalNode());


        //TODO Do some assertions after cuts if the expected costs values correspond. E.g. mutualCost[i][i]-cutValue[i], mutualCost[i][j]==currentMutualCost
        while(shiftsDone<maxMoveForEndCut&&proceedSearch){


            auto itFirst=baseEdges[pointerFirst].find(descSecond);
            bool firstExists=itFirst!=baseEdges[pointerFirst].end();
            auto itSecond=baseEdges[predFirst].find(pointerSecond);
            bool secondExists=itSecond!=baseEdges[predFirst].end();

            if(firstExists&&!secondExists){  //move the second pointer
                double cutPrice=itFirst->second+firstCutValue+secondCutShifted;
                if(cutPrice+currentMutualCost<-eps){//do the cut only if connection has improvement
                    double expectedPrimal=currentPrimalValue+getOutputCost(pointerSecond)+getInputCost(descSecond)-cutValues[pointerSecond]-baseEdges[pointerSecond].at(descSecond);
                    currentMutualCost-=lossOfSecondCost[0];
                    currentNegativeCost-=lossOfSecondCost[1];
                    double positiveMutualCost=currentMutualCost-currentNegativeCost;
                    if(positiveMutualCost<=mergeThreshold*abs(currentNegativeCost)){

                        if(pointerFirst!=lastVertices[pathIndex1]){
                            expectedPrimal+=getOutputCost(pointerFirst)+getInputCost(neighboringVertices[pointerFirst])-cutValues[pointerFirst]-baseEdges[pointerFirst].at(neighboringVertices[pointerFirst]);
                            cutOnePathWithMutualCostUpdate(pointerFirst);

                        }
                        cutOnePathWithMutualCostUpdate(pointerSecond);

                        size_t newPathIndex2=numberOfPaths-1;
                        assert(abs(cummulativeCosts[pathIndex1][newPathIndex2]-currentMutualCost)<eps);
                        assert(abs(negativeMutualCosts[pathIndex1][newPathIndex2]-currentNegativeCost)<eps);
                        baseEdgeVertices={pointerFirst,descSecond};
#ifndef DEBUG
                        currentPrimal(false);
                        assert(abs(currentPrimalValue-expectedPrimal)/abs(currentPrimalValue)<1e-10);
                        if(diagnostics()) std::cout<<"cuts done "<<std::endl;
#endif
                        //cutFound=true;
                    }
                }
                proceedSearch=false;
                //cut from first and terminate

            }
            else if(secondExists&&!firstExists){  //move the first pointer
                double cutPrice=itSecond->second+secondCutValue+firstCutShifted;
                if(cutPrice+currentMutualCost<-eps){
                    double expectedPrimal=currentPrimalValue+getInputCost(pointerFirst)+getOutputCost(predFirst)-baseEdges[predFirst].at(pointerFirst)-cutValues[predFirst];
                    currentMutualCost-=lossOfFirstCost[0];
                    currentNegativeCost-=lossOfFirstCost[1];

                    double positiveMutualCost=currentMutualCost-currentNegativeCost;
                    if(positiveMutualCost<=mergeThreshold*abs(currentNegativeCost)){

                        size_t newPathIndex2=pathIndex2;
                        if(pointerSecond!=startingVertices[pathIndex2]){

                            newPathIndex2=numberOfPaths-1;
                            expectedPrimal+=getInputCost(pointerSecond)+getOutputCost(reverseNeighbors[pointerSecond])-baseEdges[reverseNeighbors[pointerSecond]].at(pointerSecond)-cutValues[reverseNeighbors[pointerSecond]];
                            cutOnePathWithMutualCostUpdate(reverseNeighbors[pointerSecond]);
                        }
                        cutOnePathWithMutualCostUpdate(predFirst);


                        assert(abs(cummulativeCosts[pathIndex1][newPathIndex2]-currentMutualCost)<eps);
                        assert(abs(negativeMutualCosts[pathIndex1][newPathIndex2]-currentNegativeCost)<eps);
                        baseEdgeVertices={predFirst,pointerSecond};
#ifndef DEBUG
                        currentPrimal(false);
                        assert(abs(currentPrimalValue-expectedPrimal)/abs(currentPrimalValue)<1e-10);
                        if(diagnostics()) std::cout<<"cuts done "<<std::endl;
#endif
                        //cutFound=true;
                    }
                }
                proceedSearch=false;

            }
            else if(secondExists&&firstExists){
                double cutPriceShiftFirst=itSecond->second+secondCutValue+firstCutShifted;
                double cutPriceShiftSecond=itFirst->second+firstCutValue+secondCutShifted;

                if(cutPriceShiftFirst<cutPriceShiftSecond&&cutPriceShiftFirst+currentMutualCost<-eps){   //move first
                    double expectedPrimal=currentPrimalValue+getInputCost(pointerFirst)+getOutputCost(predFirst)-baseEdges[predFirst].at(pointerFirst)-cutValues[predFirst];
                    currentMutualCost-=lossOfFirstCost[0];
                    currentNegativeCost-=lossOfFirstCost[1];
                    double positiveMutualCost=currentMutualCost-currentNegativeCost;
                    if(positiveMutualCost<=mergeThreshold*abs(currentNegativeCost)){

                        size_t newPathIndex2=pathIndex2;
                        if(pointerSecond!=startingVertices[pathIndex2]){

                            newPathIndex2=numberOfPaths-1;
                            expectedPrimal+=getInputCost(pointerSecond)+getOutputCost(reverseNeighbors[pointerSecond])-baseEdges[reverseNeighbors[pointerSecond]].at(pointerSecond)-cutValues[reverseNeighbors[pointerSecond]];
                            cutOnePathWithMutualCostUpdate(reverseNeighbors[pointerSecond]);
                        }
                        cutOnePathWithMutualCostUpdate(predFirst);

                        assert(abs(cummulativeCosts[pathIndex1][newPathIndex2]-currentMutualCost)<eps);
                        assert(abs(negativeMutualCosts[pathIndex1][newPathIndex2]-currentNegativeCost)<eps);
                        baseEdgeVertices={predFirst,pointerSecond};
#ifndef DEBUG
                        currentPrimal(false);
                        assert(abs(currentPrimalValue-expectedPrimal)/abs(currentPrimalValue)<1e-10);
                        if(diagnostics()) std::cout<<"cuts done "<<std::endl;
#endif
                        // cutFound=true;
                    }

                }
                else if(cutPriceShiftSecond<cutPriceShiftFirst&&cutPriceShiftSecond+currentMutualCost<-eps){       //move second
                    double expectedPrimal=currentPrimalValue+getOutputCost(pointerSecond)+getInputCost(descSecond)-cutValues[pointerSecond]-baseEdges[pointerSecond].at(descSecond);
                    currentMutualCost-=lossOfSecondCost[0];
                    currentNegativeCost-=lossOfSecondCost[1];

                    double positiveMutualCost=currentMutualCost-currentNegativeCost;
                    if(positiveMutualCost<=mergeThreshold*abs(currentNegativeCost)){
                        if(pointerFirst!=lastVertices[pathIndex1]){

                            expectedPrimal+=getOutputCost(pointerFirst)+getInputCost(neighboringVertices[pointerFirst])-cutValues[pointerFirst]-baseEdges[pointerFirst].at(neighboringVertices[pointerFirst]);
                            cutOnePathWithMutualCostUpdate(pointerFirst);
                        }
                        cutOnePathWithMutualCostUpdate(pointerSecond);

                        size_t newPathIndex2=numberOfPaths-1;
                        assert(abs(cummulativeCosts[pathIndex1][newPathIndex2]-currentMutualCost)<eps);
                        assert(abs(negativeMutualCosts[pathIndex1][newPathIndex2]-currentNegativeCost)<eps);
                        baseEdgeVertices={pointerFirst,descSecond};
#ifndef DEBUG
                        currentPrimal(false);
                        assert(abs(currentPrimalValue-expectedPrimal)/abs(currentPrimalValue)<1e-10);
                        if(diagnostics()) std::cout<<"cuts done "<<std::endl;
#endif
                        //  cutFound=true;
                    }

                }
                proceedSearch=false;

            }
            else{
                double cutPriceShiftFirst=secondCutValue+firstCutShifted;
                double cutPriceShiftSecond=firstCutValue+secondCutShifted;


                if(cutPriceShiftFirst<cutPriceShiftSecond){

                    double extraInputCost=getInputCost(pointerFirst);
                    assert(baseEdges[predFirst].count(pointerFirst)>0);
                    double baseCost=baseEdges[predFirst].at(pointerFirst);
                    pointerFirst=predFirst;
                    firstCutValue=-cutValues[pointerFirst]+extraInputCost-baseCost;
                    currentMutualCost-=lossOfFirstCost[0];
                    currentNegativeCost-=lossOfFirstCost[1];
                    predFirst=reverseNeighbors[pointerFirst];
                    if(predFirst==pInstance->getSourceNode()){
                        break;
                    }

                }
                else{
                    double extraOutputCost=getOutputCost(pointerSecond);
                    assert(baseEdges[pointerSecond].count(descSecond)>0);
                    double baseCost=baseEdges[pointerSecond].at(descSecond);
                    secondCutValue=-cutValues[pointerSecond]+extraOutputCost-baseCost;

                    pointerSecond=descSecond;

                    currentMutualCost-=lossOfSecondCost[0];
                    currentNegativeCost-=lossOfSecondCost[1];

                    descSecond=neighboringVertices[descSecond];
                    if(descSecond==pInstance->getTerminalNode()){
                        break;
                    }


                }
                lossOfFirstCost=mutualCostToSubtract(pointerFirst,pointerSecond,pathIndex2,true); //TODO take into account what has been removed on the other side
                lossOfSecondCost=mutualCostToSubtract(pointerSecond,pointerFirst,pathIndex1,false);

                firstCutShifted=-cutValues[predFirst]-lossOfFirstCost[0]+getInputCost(pointerFirst)-baseEdges[predFirst].at(pointerFirst);
                secondCutShifted=-cutValues[pointerSecond]-lossOfSecondCost[0]+getOutputCost(pointerSecond)-baseEdges[pointerSecond].at(descSecond);
                shiftsDone++;

            }

        }
    }

    return baseEdgeVertices;

}



template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::cutToEnableConnections(){
    initLastVertices();
    computeCummulativeCostsGlobal();
    initReverseNeighbors();

    size_t noAssignment=pInstance->getTerminalNode();

    size_t originalNumberOfPaths=numberOfPaths;
    std::vector<size_t> bestNeighbors(numberOfPaths,noAssignment);
    std::vector<size_t> reverseBest(numberOfPaths,noAssignment);
    std::vector<double> bestValues(numberOfPaths,0);
    std::vector<char> validPaths(numberOfPaths,1);

    for (size_t i = 0; i < originalNumberOfPaths; ++i) {
        size_t lastVertex=lastVertices[i];
        bestNeighbors[i]=noAssignment;
        bestValues[i]=0;

        for (size_t j = 0; j < originalNumberOfPaths; ++j) {
            if(i==j) continue;
            size_t firstVertex=startingVertices[j];
            if(firstVertex!=pInstance->getTerminalNode()){
                //if(firstVertex>lastVertex){
                    auto it=baseEdges[lastVertex].find(firstVertex);
                    if(it==baseEdges[lastVertex].end()&&cummulativeCosts[i][j]<0){
                        double negativeCost=negativeMutualCosts[i][j];
                        double positiveCost=cummulativeCosts[i][j]-negativeMutualCosts[i][j];
                        //if(cummulativeCosts[i][j]<bestValues[i]&&abs(negativeCost)*mergeThreshold>=positiveCost){
                        if(cummulativeCosts[i][j]<bestValues[i]){
                              bestNeighbors[i]=j;
                              bestValues[i]=cummulativeCosts[i][j];

                        }

                    }
                    else if(it!=baseEdges[lastVertex].end()&&cummulativeCosts[i][j]<0){
                        double connectionValue=cummulativeCosts[i][j]+it->second;
                        double negativeCost=negativeMutualCosts[i][j];
                        double positiveCost=cummulativeCosts[i][j]-negativeMutualCosts[i][j];
//                        if(connectionValue<bestValues[i]&&abs(negativeCost)*mergeThreshold>=positiveCost){
                        if(connectionValue<bestValues[i]){
                            bestNeighbors[i]=j;
                            bestValues[i]=connectionValue;
                        }

                    }
               // }
            }

        }
        size_t bestNeighbor=bestNeighbors[i];
        double bestValue=bestValues[i];
        if(bestNeighbor!=noAssignment){
            if(reverseBest[bestNeighbor]!=noAssignment){
                size_t origReverse=reverseBest[bestNeighbor];
                double origValue=bestValues[origReverse];
                if(origValue>bestValue){
                    reverseBest[bestNeighbor]=i;
                    bestNeighbors[origReverse]=noAssignment;
                    bestValues[origReverse]=0;
                }
                else{
                    bestNeighbors[i]=noAssignment;
                    bestValues[i]=0;
                }
            }
            else{
                reverseBest[bestNeighbor]=i;
            }
        }

    }

    std::vector<size_t> oldToNewIndex(numberOfPaths);

    for (size_t i = 0; i < numberOfPaths; ++i) {
        oldToNewIndex[i]=i;
    }

    for (size_t i = 0; i < originalNumberOfPaths; ++i) {

        if(bestNeighbors[i]!=noAssignment){

//            size_t neighborIndex=bestNeighbors[i];
//            bestNeighbors[i]=oldToNewIndex[neighborIndex];
            size_t firstPathIndex=oldToNewIndex[i];
            size_t secondPathIndex=bestNeighbors[i];
            size_t firstVertex=startingVertices[secondPathIndex];
            size_t lastVertex=lastVertices[firstPathIndex];
            assert(firstVertex!=pInstance->getTerminalNode());
            assert(oldToNewIndex[secondPathIndex]==secondPathIndex);
            //if(firstVertex>lastVertex){
            auto it=baseEdges[lastVertex].find(firstVertex);
            if(it==baseEdges[lastVertex].end()){

                if(cummulativeCosts[firstPathIndex][secondPathIndex]<0){
                    if(diagnostics())std::cout<<"finding best cut "<<firstPathIndex<<","<<secondPathIndex<<std::endl;
                    std::array<size_t,2> baseEdge=findBestCut(firstPathIndex,secondPathIndex);
                    if(baseEdge[0]!=pInstance->getTerminalNode()){
                        size_t newLabel=vertexLabels[baseEdge[1]];
                        assert(newLabel==numberOfPaths||newLabel==secondPathIndex+1);
                        oldToNewIndex[secondPathIndex]=newLabel-1;
                        auto itBE=baseEdges[baseEdge[0]].find(baseEdge[1]);
                        assert(itBE!=baseEdges[baseEdge[0]].end());
                    }
                }

            }
        }
    }


}



template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::mergePaths(){


    bool tryNext=true;

    while(tryNext){
        size_t bestStartingPathIndex=numberOfPaths; //Probably using multimap <double,size_t>, always update all paths that need to be updated? Or for over paths and find best during update
        double bestStartingPathValue=0;
        std::vector<size_t> bestNeighboringPaths(numberOfPaths,numberOfPaths);

        double missedConnectionsCost=0;
        size_t missedConnections=0;


        std::vector<std::map<size_t,double>> costsBetweenPaths(numberOfPaths);
        assert(numberOfPaths==lastVertices.size());
        for (size_t i = 0; i < numberOfPaths; ++i) {
            size_t lastVertex=lastVertices[i];
            size_t bestNeighbor=numberOfPaths;
            double bestNeighborValue=0;
            for (size_t j = 0; j < numberOfPaths; ++j) {
                size_t firstVertex=startingVertices[j];
                if(firstVertex!=pInstance->getTerminalNode()){
                    if(firstVertex>lastVertex){
                        auto it=baseEdges[lastVertex].find(firstVertex);
                        if(it!=baseEdges[lastVertex].end()){
                            double cost=it->second;
                            cost-=getInputCost(j);
                            cost-=getOutputCost(i);
                            cost+=cummulativeCosts[i][j];
                            costsBetweenPaths[i][j]=cost;

                            if(cost<bestNeighborValue){

                                double negativeCost=negativeMutualCosts[i][j];
                                double positiveCost=cummulativeCosts[i][j]-negativeMutualCosts[i][j];
                                if(abs(negativeCost)*mergeThreshold>=positiveCost){

                                    bestNeighborValue=cost;
                                    bestNeighbor=j;
                                }

                            }
                        }

                    }
                }

            }
            bestNeighboringPaths[i]=bestNeighbor;
            if(bestNeighborValue<bestStartingPathValue){
                bestStartingPathValue=bestNeighborValue;
                bestStartingPathIndex=i;
            }
        }


        if(bestStartingPathIndex<numberOfPaths){
            assert(bestNeighboringPaths[bestStartingPathIndex]<numberOfPaths);
            connectToOnePath(bestStartingPathIndex, costsBetweenPaths,bestNeighboringPaths[bestStartingPathIndex]);

            double primal=currentPrimal(false);
            if(diagnostics()) std::cout<<"primal after merge "<<primal<<std::endl;

        }
        else{
            tryNext=false;
        }



    }


}

template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::initLastVertices(){
    lastVertices=std::vector<size_t>(startingVertices.size());
    for (size_t i = 0; i < neighboringVertices.size(); ++i) {
        if(neighboringVertices[i]==pInstance->getTerminalNode()&&vertexLabels[i]!=0){
            size_t pathIndex=vertexLabels[i]-1;
            assert(lastVertices[pathIndex]==0);
            lastVertices[pathIndex]=i;
        }

    }

}


template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::initCutValues(){  //get lifted values of cutting behind a vertex
    cutValues=std::vector<double> (pInstance->getNumberOfVertices()-2);

    const LdpDirectedGraph& liftedGraph=getLiftedGraph();


    //Addding new lifted edges to cummulative costs
    for (int i = 0; i < liftedGraph.getNumberOfVertices(); ++i) {
        size_t label=vertexLabels[i];
        if(label>0){
            auto iter=liftedGraph.forwardNeighborsBegin(i);
            for (;iter!=liftedGraph.forwardNeighborsEnd(i);iter++) {
                size_t vertex2=iter->first;
                size_t label2=vertexLabels[vertex2];
                if(label2==label){

                    double cost=iter->second;
                    size_t vertex1=i;
                    while(vertex1!=vertex2){
                        cutValues[vertex1]+=cost;
                        vertex1=neighboringVertices[vertex1];
                        assert(vertex1!=pInstance->getTerminalNode());
                    }
                }
            }
        }
    }


}



template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::updateCutValues(const std::vector<size_t>& oldLabels, const std::vector<size_t>& oldNeighbors){  //get lifted values of cutting behind a vertex
   // std::vector<double> cutValues(pInstance->getNumberOfVertices()-2);

    const LdpDirectedGraph& liftedGraph=getLiftedGraph();


    for (int i = 0; i < liftedGraph.getNumberOfVertices(); ++i) {
        size_t label=oldLabels[i];
        size_t newLabel=vertexLabels[i];
        if(label>0){
            auto iter=liftedGraph.forwardNeighborsBegin(i);
            for (;iter!=liftedGraph.forwardNeighborsEnd(i);iter++) {
                size_t vertex2=iter->first;
                size_t label2=oldLabels[vertex2];
                size_t newLabel2=vertexLabels[vertex2];
                if(label2==label&&newLabel!=newLabel2){
                    double cost=iter->second;
                    size_t vertex1=i;
                    while(vertex1!=vertex2){
                        cutValues[vertex1]-=cost;

                        vertex1=oldNeighbors[vertex1];
                        assert(vertex1!=pInstance->getTerminalNode());
                    }
                }
            }
        }
    }


}



template<class SNC_FACTOR>
std::vector<size_t> LdpPrimalHeuristics<SNC_FACTOR>::findCutCandidates(const std::vector<double>& cutValues){
    double minValueToCut=0;
    //std::vector<double> cutValues=initCandiadteCuts();
    std::vector<size_t> bestCutsInPaths(numberOfPaths,pInstance->getTerminalNode());
    for (size_t i = 0; i < numberOfPaths; ++i) {
        size_t currentVertex=startingVertices[i];
        if(currentVertex!=pInstance->getTerminalNode()){
            double bestCutValue=minValueToCut;
            while(currentVertex!=pInstance->getTerminalNode()){
                double cutValue=cutValues[currentVertex];
                size_t nextVertex=neighboringVertices[currentVertex];
                assert(baseEdges[currentVertex].count(nextVertex)>0);
                if(nextVertex!=pInstance->getTerminalNode()){
                    cutValue+=baseEdges[currentVertex][nextVertex]-getInputCost(nextVertex)-getOutputCost(currentVertex);

                    if(cutValue>bestCutValue&&nextVertex!=pInstance->getTerminalNode()){
                        bestCutsInPaths[i]=currentVertex;
                        bestCutValue=cutValue;
                    }
                }
                currentVertex=nextVertex;
            }
        }

    }
    return bestCutsInPaths;

}




template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::cutPaths(const std::vector<size_t>& cutCandidates,const std::vector<size_t> indicesOfPathsToCut){
    size_t newNumberOfPaths=numberOfPaths;
    std::vector<size_t> newLabels=vertexLabels;
    for (size_t p = 0; p < indicesOfPathsToCut.size(); ++p) {
        size_t i=indicesOfPathsToCut[p];
        assert(cutCandidates[i]!=pInstance->getTerminalNode());
        size_t cutVertex=cutCandidates[i];
        size_t nextVertex=neighboringVertices[cutVertex];
        assert(nextVertex!=pInstance->getTerminalNode());
        assert(vertexLabels[cutVertex]==i+1);
        assert(vertexLabels[nextVertex]==i+1);
        neighboringVertices[cutVertex]=pInstance->getTerminalNode();
        startingVertices.push_back(nextVertex);
        size_t newLabel=startingVertices.size();
        while(nextVertex!=pInstance->getTerminalNode()){
            newLabels[nextVertex]=newLabel;
            nextVertex=neighboringVertices[nextVertex];
        }

    }
    numberOfPaths=startingVertices.size();
    vertexLabels=newLabels;


}



template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::improveAllPaths(){

    double primal=currentPrimal();
    //std::cout<<"primal value "<<primal<<std::endl;

    bool tryCut=true;
    initCutValues();
    while(tryCut){
        std::vector<size_t> cutCandidates=findCutCandidates(cutValues);
        std::vector<size_t> indicesOfPathsToCut;
        double improvement=0;
        for (size_t i = 0; i < numberOfPaths; ++i) {

            if(cutCandidates[i]!=pInstance->getTerminalNode()){
                indicesOfPathsToCut.push_back(i);
                size_t cutVertex=cutCandidates[i];
                double value=cutValues[cutVertex];
                size_t nextVertex=neighboringVertices[cutVertex];
                value+=baseEdges[cutVertex][nextVertex]-getInputCost(nextVertex)-getOutputCost(cutVertex);
                improvement+=value;
                if(debug())std::cout<<"improvement from path "<<i<<": "<<value<<", cut vertex: "<<cutVertex<<std::endl;

            }
        }
        double expectedPrimal=primal-improvement;
        if(diagnostics()) std::cout<<"expected new primal "<<(primal-improvement)<<std::endl;
        if(indicesOfPathsToCut.size()>0){
            std::vector<size_t> oldVertexLabels=vertexLabels;
            std::vector<size_t> oldNeighbors=neighboringVertices;
            cutPaths(cutCandidates,indicesOfPathsToCut);
            primal=currentPrimal();
            assert(std::abs(primal-expectedPrimal)/abs(expectedPrimal)<1e-10);
            if(diagnostics()) std::cout<<"primal value "<<primal<<std::endl;
            updateCutValues(oldVertexLabels,oldNeighbors);

        }
        else{
            tryCut=false;
        }
    }


    cutToEnableConnections();
    if(diagnostics()) std::cout<<"call merge paths"<<std::endl;
    mergePaths();

    finalizeResults(true);


}




template<class SNC_FACTOR>
double LdpPrimalHeuristics<SNC_FACTOR>::currentPrimal(bool checkImprovement){
    const LdpDirectedGraph& liftedGraph=getLiftedGraph();
    const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
    costsOfPaths=std::vector<double>(numberOfPaths);

    size_t counterPaths=0;
    double newPrimalValue=0;

    assert(startingVertices.size()==numberOfPaths);
    for (int i = 0; i < numberOfPaths; ++i) {
        double pathCost=0;
        if(startingVertices[i]!=pInstance->getTerminalNode()){

            counterPaths++;
            size_t currentVertex=startingVertices[i];
            pathCost+=getInputCost(currentVertex);   //output cost is included in the while cycle
            newPrimalValue+=getInputCost(currentVertex);

            std::vector<size_t> newPath;

            while(currentVertex!=pInstance->getTerminalNode()){
               // vertexLabels[currentVertex]=counterPaths;
                newPath.push_back(currentVertex);
                size_t newVertex=neighboringVertices[currentVertex];
                double beCost=baseEdges[currentVertex].at(newVertex);
                newPrimalValue+=beCost;
                pathCost+=beCost;
                if(newVertex!=pInstance->getTerminalNode()){
                    assert(vertexLabels[newVertex]==i+1);
                }
                currentVertex=newVertex;

            }
           // adjustedPaths.push_back(newPath);
        }
        costsOfPaths[i]=pathCost;
    }


    std::vector<double> costOfPathsBase=costsOfPaths;
    std::vector<double> costOfPathsLifted(numberOfPaths);

    double edgeLoss=0;
    for (int i = 0; i < liftedGraph.getNumberOfVertices(); ++i) {
        if(vertexLabels[i]!=0){
            for (auto iter=liftedGraph.forwardNeighborsBegin(i);iter!=liftedGraph.forwardNeighborsEnd(i);iter++) {
                size_t vertex2=iter->first;

                if(vertexLabels[vertex2]==vertexLabels[i]){
                    costsOfPaths[vertexLabels[vertex2]-1]+=iter->second;
                    costOfPathsLifted[vertexLabels[vertex2]-1]+=iter->second;
                    newPrimalValue+=iter->second;

                }
            }
        }
    }

#ifndef NDEBUG
    double checkValue=0;
    for (size_t i = 0; i < numberOfPaths; ++i) {
        checkValue+=costsOfPaths[i];
    }
    assert(abs(checkValue-newPrimalValue)/abs(newPrimalValue)<1e-10);

#endif

    assert(!checkImprovement||newPrimalValue<=currentPrimalValue+eps);
    currentPrimalValue=newPrimalValue;
    return newPrimalValue;
   // std::cout<<"current primal value "<<newPrimalValue<<std::endl;


}




template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::finalizeResults(bool changeSNC){
    currentPrimal(true);
    const LdpDirectedGraph& liftedGraph=getLiftedGraph();
    const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
    std::vector<std::array<SNC_FACTOR*,2>>& single_node_cut_factors_=*pSNCFactors;
    if(changeSNC){

          for (std::size_t graph_node = 0; graph_node < liftedGraph.getNumberOfVertices(); ++graph_node) {
              single_node_cut_factors_[graph_node][0]->get_factor()->setNoBaseEdgeActive();
              single_node_cut_factors_[graph_node][1]->get_factor()->setNoBaseEdgeActive();
          }
    }


    size_t counterPaths=0;
    double newPrimalValue=0;
    std::vector<size_t> newStartingVertices;
    vertexLabels=std::vector<size_t>(pInstance->getNumberOfVertices()-2);
    assert(startingVertices.size()==numberOfPaths);
    for (int i = 0; i < numberOfPaths; ++i) {
        if(startingVertices[i]!=pInstance->getTerminalNode()&&costsOfPaths[i]<=0){
        //if(startingVertices[i]!=pInstance->getTerminalNode()){
            //newPrimalValue+=pInstance->parameters.getOutputCost();
            newPrimalValue+=getInputCost(startingVertices[i]);
            //newPrimalValue+=pInstance->parameters.getInputCost();
            counterPaths++;
            size_t currentVertex=startingVertices[i];
            //newStartingVertices.push_back(currentVertex);
            std::vector<size_t> newPath;
            if(changeSNC){
                auto* pSNC=single_node_cut_factors_[currentVertex][0]->get_factor();
                pSNC->setBaseEdgeActiveWithID(pInstance->getSourceNode());
            }
            while(currentVertex!=pInstance->getTerminalNode()){
                vertexLabels[currentVertex]=counterPaths;
                newPath.push_back(currentVertex);
                size_t newVertex=neighboringVertices[currentVertex];
                if(changeSNC){
                    auto* pSNCOut=single_node_cut_factors_[currentVertex][1]->get_factor();
                    pSNCOut->setBaseEdgeActiveWithID(newVertex);

                }
                if(newVertex!=pInstance->getTerminalNode()){
                    if(changeSNC){
                        auto* pSNCIn=single_node_cut_factors_[newVertex][0]->get_factor();
                        pSNCIn->setBaseEdgeActiveWithID(currentVertex);
                    }
                }
                newPrimalValue+=baseEdges[currentVertex].at(newVertex);
                currentVertex=newVertex;

            }
            adjustedPaths.push_back(newPath);
            newStartingVertices.push_back(startingVertices[i]);
        }

    }
    startingVertices=newStartingVertices;
    assert(counterPaths==startingVertices.size());
    assert(adjustedPaths.size()==counterPaths);
    numberOfPaths=counterPaths;

    for (int i = 0; i < liftedGraph.getNumberOfVertices(); ++i) {
        if(vertexLabels[i]!=0){
            for (auto iter=liftedGraph.forwardNeighborsBegin(i);iter!=liftedGraph.forwardNeighborsEnd(i);iter++) {
                size_t vertex2=iter->first;
                if(vertexLabels[vertex2]==vertexLabels[i]){
                    newPrimalValue+=iter->second;
                }
            }
        }
    }

    currentPrimalValue=newPrimalValue;
    if (diagnostics()) std::cout<<"new primal value "<<newPrimalValue<<std::endl;


}




}
#endif // LDP_PRIMAL_HEURISTICS_HXX
