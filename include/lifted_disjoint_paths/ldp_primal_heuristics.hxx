#ifndef LDP_PRIMAL_HEURISTICS_HXX
#define LDP_PRIMAL_HEURISTICS_HXX

#include <cstdlib>
#include <vector>
#include <lifted_disjoint_paths/ldp_instance.hxx>
#include "graph_matching/graph_matching_input.h"
#include "graph_matching/min_cost_flow_factor_ssp.hxx"
#include "graph_matching/matching_problem_input.h"

namespace LPMP{

template <class SNC_FACTOR>
class LdpPrimalHeuristics{
public:
    LdpPrimalHeuristics(const std::vector<size_t>& _vertexLabels,
                        const std::vector<size_t>& _startingVertices,
                        //const std::vector<std::vector<size_t>>& _paths,
                        const std::vector<size_t>& _neighboringVertices,
                        const lifted_disjoint_paths::LdpInstance* _pInstance,
                        std::vector<std::array<SNC_FACTOR*,2>>* _p_single_node_cut_factors)
      {
        vertexLabels=_vertexLabels;
        startingVertices=_startingVertices;
       // paths=_paths;
        neighboringVertices=_neighboringVertices;
        pInstance=_pInstance;
        pSNCFactors=_p_single_node_cut_factors;

        numberOfPaths=startingVertices.size();
        const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();

        baseEdges=std::vector<std::map<size_t,double>>(baseGraph.getNumberOfVertices());
        for (int i = 0; i < baseGraph.getNumberOfVertices(); ++i) {
            for (auto iter=baseGraph.forwardNeighborsBegin(i);iter!=baseGraph.forwardNeighborsEnd(i);iter++) {
                baseEdges[i][iter->first]=iter->second;
            }

        }

        double pValue=currentPrimal();
        if(diagnostics()) std::cout<<"init primal value "<<pValue<<std::endl;

       // cutAllPaths();
        pValue=currentPrimal();
        if(diagnostics()) std::cout<<"init primal value "<<pValue<<std::endl;


//        cummulativeCosts=std::vector<std::vector<double>>(numberOfPaths);
//        for (int i = 0; i < numberOfPaths; ++i) {
//            cummulativeCosts[i]=std::vector<double>(numberOfPaths);
//        }

        //pointersToLiftedNeighbors=std::vector<size_t>(pInstance->getNumberOfVertices()-2);

        //pointersToPaths=std::vector<size_t>(numberOfPaths,-1); //some paths have not started yet
        currentTime=1;

        //firstRelevantStartVertex=0;

       pointersToPaths=startingVertices;

       currentPrimalValue=0;

    }

    ~LdpPrimalHeuristics(){

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

void cutAllPaths();


private:
    void finalizeResults(bool changeSNC);
    void evaluateTimeCuts();
    double currentPrimal();

    void prepareCummulativeCostsInCurrentTime();

    void initLastVertices();
    void mergePaths();
    void connectTwoPaths(size_t indexFirst,size_t indexLast);
    void connectToOnePath(size_t indexOfPath, const std::vector<std::map<size_t,double>>& costsBetweenPaths, size_t bestNeighborIndex);

    void computeCummulativeCostsGlobal();
    std::vector<double> initCutValues();
    void updateCutValues(const std::vector<size_t>& oldLabels, const std::vector<size_t>& oldNeighbors, std::vector<double>& cutValues);
    std::vector<size_t> findCutCandidates(const std::vector<double>& cutValues);
    void cutPaths(const std::vector<size_t>& cutCandidates,const std::vector<size_t> indicesOfPathsToCut);

    void initCummulativeCosts();
    std::vector<std::array<int,2>> setPathPointers();
    double initLAP(LPMP::linear_assignment_problem_input& lapInput, std::vector<std::array<int, 2> > &baseIndices);

    double getBaseEdgeCost(size_t v,size_t w)const{
         const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
        const std::vector<std::array<SNC_FACTOR*,2>>& single_node_cut_factors_=*pSNCFactors;
        auto pOutFactor=single_node_cut_factors_[v][1]->get_factor();
        const std::unordered_map<size_t,size_t>& baseIDsToIndices=pOutFactor->getBaseIDsToIndices();

        auto it=baseIDsToIndices.find(w);
        assert(it!=baseIDsToIndices.end());
        size_t index=it->second;
        double cost=baseGraph.getForwardEdgeCost(v,index);
        return cost;


     }



    std::vector<std::map<size_t,double>> baseEdges;
    std::vector<size_t> vertexLabels;
    std::vector<size_t> startingVertices;
    std::vector<std::vector<size_t>> adjustedPaths;
    std::vector<size_t> neighboringVertices;
    const lifted_disjoint_paths::LdpInstance* pInstance;
    size_t numberOfPaths;
    double currentPrimalValue;

    std::vector<size_t> lastVertices;
   // const std::vector<std::array<SNC_FACTOR*,2>>& single_node_cut_factors_;
    std::vector<std::array<SNC_FACTOR*,2>>* pSNCFactors;

    size_t currentTime;
   // size_t firstRelevantStartVertex;
    //size_t lastRelevantStartVertex; //plus max time gap time frames from first start vertex
    std::vector<size_t> pointersToLiftedNeighbors; //size of number of vertices, pointer to lifted neighbors relevant for current inspected time
    //pointer (vertex index) to last vertex in paths relevant to current timeFrame

     std::vector<size_t> pointersToPaths;
    std::vector<std::vector<double>> cummulativeCosts;  //dimension number of paths x number of paths, really the full dimension!

};


template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::computeCummulativeCostsGlobal(){
    initCummulativeCosts();

    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();


    //Addding new lifted edges to cummulative costs
    for (int i = 0; i < liftedGraph.getNumberOfVertices(); ++i) {
        size_t label=vertexLabels[i];
        if(label>0){
            auto iter=liftedGraph.forwardNeighborsBegin(i);
            for (;iter!=liftedGraph.forwardNeighborsEnd(i);iter++) {
                size_t vertex2=iter->first;
                size_t label2=vertexLabels[vertex2];
//                if(label==1&&label2==51){
//                    std::cout<<"from 0 to 50, edge "<<i<<","<<vertex2<<":"<<iter->second<<std::endl;
                //                }
                if(label2>0){
                    assert(label2>0);
                    assert(label-1<cummulativeCosts.size());
                    assert(label2-1<cummulativeCosts[label-1].size());
                    cummulativeCosts[label-1][label2-1]+=iter->second;
                }
            }
        }
    }

}



template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::connectToOnePath(size_t indexOfPath, const std::vector<std::map<size_t,double>>& costsBetweenPaths,size_t bestNeighborIndex){

    const VertexGroups<>& vg=(pInstance->vertexGroups);

    size_t maxTimeToUse=vg.getGroupIndex(startingVertices[bestNeighborIndex]);

    //The following for can be already part of DFS
    std::stack<size_t> stackOfPathIndices;
    stackOfPathIndices.push(indexOfPath);
    std::stack<double> stackOfCosts;
    std::stack<size_t> stackOfPredecessors;
    stackOfCosts.push(0.0);
    std::vector<size_t> currentPath;
    std::vector<size_t> bestPathSoFar={indexOfPath};
    double bestCostSoFar=costsBetweenPaths[indexOfPath].at(bestNeighborIndex);


    while(!stackOfPathIndices.empty()){
        //std::vector<size_t> alternativeNeighbors;
        size_t currentVertex=stackOfPathIndices.top();
        double currentPathCost=stackOfCosts.top();
        const std::map<size_t,double>& neighbors=costsBetweenPaths[currentVertex];
        currentPath.push_back(currentVertex);

        if(!stackOfPredecessors.empty()&&stackOfPredecessors.top()==currentVertex){
            auto it=neighbors.find(bestNeighborIndex);
            if(it!=neighbors.end()){
                double costToLast=it->second;
                if(costToLast<0){
                    costToLast+=currentPathCost;
                    if(costToLast<bestCostSoFar){
                        bestPathSoFar=currentPath;
                    }
                }
            }
            stackOfPredecessors.pop();
            stackOfCosts.pop();
            stackOfPathIndices.pop();
            currentPath.pop_back();

        }
        else{
            stackOfPredecessors.push(currentVertex);
            for (auto iter=neighbors.begin();iter!=neighbors.end();iter++) {
                if(iter->second<0&&iter->first!=bestNeighborIndex){
                    size_t middleNeighbor=iter->first;
                    auto it=costsBetweenPaths[middleNeighbor].find(bestNeighborIndex);
                    if(it!=costsBetweenPaths[middleNeighbor].end()&&it->second<0){
                        stackOfPathIndices.push(middleNeighbor);
                        stackOfCosts.push(currentPathCost+iter->second);
                    }
                    else{
                        size_t timeMiddle=vg.getGroupIndex(lastVertices[middleNeighbor]);
                        if(timeMiddle+1<maxTimeToUse&&cummulativeCosts[middleNeighbor][bestNeighborIndex]<0){
                            stackOfPathIndices.push(middleNeighbor);
                            stackOfCosts.push(currentPathCost+iter->second);
                        }
                    }
                    //maybe add anoter option whre middle neighbor and best neighbor are not connected but reachable?
                }
            }
        }


    }
    bestPathSoFar.push_back(bestNeighborIndex);

    if(diagnostics()){
        std::cout<<"merging path "<<indexOfPath<<" with "<<bestNeighborIndex<<"usning "<<std::endl;
        for (int i = 0; i < bestPathSoFar.size(); ++i) {
            std::cout<<bestPathSoFar[i]<<",";
        }
        std::cout<<std::endl;
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

    for (size_t j = 1; j < bestPathSoFar.size(); ++j) {
        size_t oldPathIndex=bestPathSoFar[j];
        std::fill(cummulativeCosts[oldPathIndex].begin(),cummulativeCosts[oldPathIndex].end(),0);

    }

    for (size_t i = 0; i < numberOfPaths; ++i) {
        for (size_t j = 1; j < bestPathSoFar.size(); ++j) {
            size_t oldPathIndex=bestPathSoFar[j];
            cummulativeCosts[i][indexOfPath]+=cummulativeCosts[i][oldPathIndex];
            cummulativeCosts[i][oldPathIndex]=0;

        }

    }



}


template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::mergePaths(){
    initLastVertices();
    computeCummulativeCostsGlobal();


    //  std::set<std::tuple<double,size_t,size_t>> queueOfEdges;
    //TODO remember best neighbor of each last vertex, plus set of best starting vertex candidates
    //The first is the costs between paths, second map <cost,path index>. Or the second always find.

    bool tryNext=true;

    while(tryNext){
        size_t bestStartingPathIndex=numberOfPaths; //Probably using multimap <double,size_t>, always update all paths that need to be updated? Or for over paths and find best during update
        double bestStartingPathValue=0;
        std::vector<size_t> bestNeighboringPaths(numberOfPaths,numberOfPaths);


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
                            cost+=cummulativeCosts[i][j];
                            costsBetweenPaths[i][j]=cost;

                            if(cost<bestNeighborValue){
                                bestNeighborValue=cost;
                                bestNeighbor=j;

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

            double primal=currentPrimal();
            if(diagnostics()) std::cout<<"primal after merge "<<primal<<std::endl;
            //TODO merge paths, update cummulative costs, best neigbor, cost between paths etc., labels, starting last vertices, neighbors
        }
        else{
            tryNext=false;
        }



    }

    //TODO after implementing path merge? change to greedy with shorter edge check!
    //Maybe: order good neighboring paths according to starting vertices, try to find possible connections
    //DFS: try nearest, from nerest nearest...
    //Greedy can select best starting path, then dfs on its neighbors (but only withint the negative neighbors), other connections in next iterations







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
std::vector<double> LdpPrimalHeuristics<SNC_FACTOR>::initCutValues(){  //get lifted values of cutting behind a vertex
    std::vector<double> cutValues(pInstance->getNumberOfVertices()-2);

    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();


    //Addding new lifted edges to cummulative costs
    for (int i = 0; i < liftedGraph.getNumberOfVertices(); ++i) {
        size_t label=vertexLabels[i];
        if(label>0){
            auto iter=liftedGraph.forwardNeighborsBegin(i);
            for (;iter!=liftedGraph.forwardNeighborsEnd(i);iter++) {
                size_t vertex2=iter->first;
                size_t label2=vertexLabels[vertex2];
                if(label2==label){
//                    if(label2==1){
//                        std::cout<<"within path 0: "<<i<<","<<vertex2<<": "<<iter->second<<std::endl;
//                    }
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
    return cutValues;

}



template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::updateCutValues(const std::vector<size_t>& oldLabels, const std::vector<size_t>& oldNeighbors, std::vector<double>& cutValues){  //get lifted values of cutting behind a vertex
   // std::vector<double> cutValues(pInstance->getNumberOfVertices()-2);

    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();


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
//                        if(vertex1==975){
//                            std::cout<<"subtract from 975 "<<cost<<", cut value "<<cutValues[vertex1]<<std::endl;

//                        }

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
                cutValue+=baseEdges[currentVertex][nextVertex];

                if(cutValue>bestCutValue&&nextVertex!=pInstance->getTerminalNode()){
                    bestCutsInPaths[i]=currentVertex;
                    bestCutValue=cutValue;
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
void LdpPrimalHeuristics<SNC_FACTOR>::cutAllPaths(){

    double primal=currentPrimal();
    //std::cout<<"primal value "<<primal<<std::endl;

    bool tryCut=true;
    std::vector<double> cutValues=initCutValues();
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
                value+=baseEdges[cutVertex][nextVertex];
                improvement+=value;
                std::cout<<"improvement from path "<<i<<": "<<value<<", cut vertex: "<<cutVertex<<std::endl;

            }
        }
        if(diagnostics()) std::cout<<"expected new primal "<<(primal-improvement)<<std::endl;
        if(indicesOfPathsToCut.size()>0){
            std::vector<size_t> oldVertexLabels=vertexLabels;
            std::vector<size_t> oldNeighbors=neighboringVertices;
            cutPaths(cutCandidates,indicesOfPathsToCut);
            primal=currentPrimal();
            if(diagnostics()) std::cout<<"primal value "<<primal<<std::endl;
            updateCutValues(oldVertexLabels,oldNeighbors,cutValues);

        }
        else{
            tryCut=false;
        }
    }

    if(diagnostics()) std::cout<<"call merge paths"<<std::endl;
    mergePaths();

    finalizeResults(true);


}




template<class SNC_FACTOR>
double LdpPrimalHeuristics<SNC_FACTOR>::currentPrimal(){
    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();
    const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
    std::vector<double> costsOfPaths(numberOfPaths);

    size_t counterPaths=0;
    double newPrimalValue=0;
    std::vector<size_t> newStartingVertices;
    assert(startingVertices.size()==numberOfPaths);
    for (int i = 0; i < numberOfPaths; ++i) {
        double pathCost=0;
        if(startingVertices[i]!=pInstance->getTerminalNode()){
            newPrimalValue+=pInstance->parameters.getOutputCost();
            newPrimalValue+=pInstance->parameters.getInputCost();
            pathCost+=pInstance->parameters.getOutputCost();
            pathCost+=pInstance->parameters.getInputCost();
            counterPaths++;
            size_t currentVertex=startingVertices[i];
            newStartingVertices.push_back(currentVertex);
            std::vector<size_t> newPath;

            while(currentVertex!=pInstance->getTerminalNode()){
               // vertexLabels[currentVertex]=counterPaths;
                newPath.push_back(currentVertex);
                size_t newVertex=neighboringVertices[currentVertex];
                if(newVertex!=pInstance->getTerminalNode()){
                    double beCost=getBaseEdgeCost(currentVertex,newVertex);
                    newPrimalValue+=beCost;
                    pathCost+=beCost;
                }
                currentVertex=newVertex;

            }
           // adjustedPaths.push_back(newPath);
        }
        costsOfPaths[i]=pathCost;
    }
    //startingVertices=newStartingVertices;
    //assert(counterPaths==startingVertices.size());
    //assert(adjustedPaths.size()==counterPaths);
    //numberOfPaths=counterPaths;

    double edgeLoss=0;
    for (int i = 0; i < liftedGraph.getNumberOfVertices(); ++i) {
        if(vertexLabels[i]!=0){
            for (auto iter=liftedGraph.forwardNeighborsBegin(i);iter!=liftedGraph.forwardNeighborsEnd(i);iter++) {
                size_t vertex2=iter->first;
//                if(vertexLabels[i]==45&&vertexLabels[vertex2]==52){
//                     std::cout<<"lost edge "<<i<<" "<<vertex2<<": "<<iter->second<<std::endl;
//                     edgeLoss+=iter->second;
//                }
                if(vertexLabels[vertex2]==vertexLabels[i]){
                    costsOfPaths[vertexLabels[vertex2]-1]+=iter->second;
                    newPrimalValue+=iter->second;
//                    if(vertexLabels[vertex2]==52){
//                       // std::cout<<"edge "<<i<<" "<<vertex2<<": "<<iter->second<<std::endl;
//                    }
                }
            }
        }
    }

//    for (int i = 0; i < numberOfPaths; ++i) {
//        std::cout<<"cost of path "<<i<<":"<<costsOfPaths[i]<<std::endl;

//    }

    assert(newPrimalValue<=currentPrimalValue+eps);
    currentPrimalValue=newPrimalValue;
    return newPrimalValue;
   // std::cout<<"current primal value "<<newPrimalValue<<std::endl;


}




template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::finalizeResults(bool changeSNC){
    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();
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
    assert(startingVertices.size()==numberOfPaths);
    for (int i = 0; i < numberOfPaths; ++i) {
        if(startingVertices[i]!=pInstance->getTerminalNode()){
            newPrimalValue+=pInstance->parameters.getOutputCost();
            newPrimalValue+=pInstance->parameters.getInputCost();
            counterPaths++;
            size_t currentVertex=startingVertices[i];
            newStartingVertices.push_back(currentVertex);
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
                    newPrimalValue+=getBaseEdgeCost(currentVertex,newVertex);
                }
                currentVertex=newVertex;

            }
            adjustedPaths.push_back(newPath);
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



template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::evaluateAll(){
    size_t iterations=pInstance->parameters.getPrimalHeuristicIterations();
    double previousValue=0;
    for (int j = 0; j < iterations; ++j) {

        for (int i = 1; i < pInstance->vertexGroups.getMaxTime(); ++i) {
            currentTime=i;
            evaluateTimeCuts();
        }
        pointersToPaths=startingVertices;
        double primalValue=currentPrimal();
        if(abs(primalValue-previousValue)/abs(primalValue)<1e-10){
            break;
        }

        previousValue=primalValue;
        if(diagnostics()) std::cout<<"primal value in iteration "<<j<<": "<<primalValue<<std::endl;
    }
    finalizeResults(true);

}

template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::initCummulativeCosts(){
    cummulativeCosts=std::vector<std::vector<double>>(numberOfPaths);
    for (int i = 0; i < numberOfPaths; ++i) {
        cummulativeCosts[i]=std::vector<double>(numberOfPaths);
    }
}

template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::prepareCummulativeCostsInCurrentTime(){


    const VertexGroups<>& vg=pInstance->getVertexGroups();
    const std::vector<size_t>& verticesInTimeLayer=vg.getGroupVertices(currentTime);
    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();


    //Addding new lifted edges to cummulative costs
    for (int i = 0; i < verticesInTimeLayer.size(); ++i) {
        size_t vertex=verticesInTimeLayer[i];
        size_t label=vertexLabels[vertex];
        if(label>0){
            auto iter=liftedGraph.forwardNeighborsBegin(vertex);
            for (;iter!=liftedGraph.forwardNeighborsEnd(vertex);iter++) {
                size_t vertex2=iter->first;
                size_t label2=vertexLabels[vertex2];
                if(label2!=0){
                    double cost=iter->second;
                    cummulativeCosts[label-1][label2-1]+=cost;
                }

            }
            auto iterBack=liftedGraph.backwardNeighborsBegin(vertex);
            for (;iterBack!=liftedGraph.backwardNeighborsEnd(vertex);iterBack++) {
                size_t vertex2=iterBack->first;
                size_t label2=vertexLabels[vertex2];
                if(label2!=0){
                    double cost=iterBack->second;
                    cummulativeCosts[label2-1][label-1]-=cost;
                }
            }
        }
    }

}


template<class SNC_FACTOR>
std::vector<std::array<int,2>> LdpPrimalHeuristics<SNC_FACTOR>::setPathPointers(){
    const VertexGroups<>& vg=pInstance->getVertexGroups();
    const std::vector<size_t>& verticesInTimeLayer=vg.getGroupVertices(currentTime);

       std::vector<std::array<int,2>> baseIndices(numberOfPaths);
    //maybe keep the option to connect in any time layer with already terminated paths or with paths that have not started yet
    for (size_t i = 0; i < pointersToPaths.size(); ++i) {
        size_t vertex=pointersToPaths[i];
        baseIndices[i]={-1,-1};
        if(vg.getGroupIndex(vertex)<=currentTime){
            size_t nextVertex=neighboringVertices[vertex];
            assert(vg.getGroupIndex(nextVertex)>=currentTime);
            if(vg.getGroupIndex(nextVertex)==currentTime){
                vertex=nextVertex;
            }
            pointersToPaths[i]=vertex;
            baseIndices[i][0]=vertex;
            nextVertex=neighboringVertices[vertex];
            if(nextVertex!=pInstance->getTerminalNode()){
                assert(vg.getGroupIndex(nextVertex)>currentTime);
                baseIndices[i][1]=nextVertex;

            }

        }
        else{ //for paths that have not started yet
            if(vertex!=pInstance->getTerminalNode()){
                baseIndices[i][1]=vertex;
            }
        }

    }

    return baseIndices;

}


template<class SNC_FACTOR>
double LdpPrimalHeuristics<SNC_FACTOR>::initLAP(LPMP::linear_assignment_problem_input& lapInput,std::vector<std::array<int,2>>& baseIndices){

    const std::vector<std::array<SNC_FACTOR*,2>>& single_node_cut_factors_=*pSNCFactors;
    const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();
    lapInput.no_left_nodes=numberOfPaths;
    lapInput.no_right_nodes=numberOfPaths;

    std::vector<std::vector<double>> pairs(numberOfPaths+1);
    for (int i = 0; i < numberOfPaths+1; ++i) {
        pairs[i]=std::vector<double>(numberOfPaths+1);
    }

    double origValue=0;

    for (int i = 0; i < numberOfPaths; ++i) {
         //lapInput.add_assignment(i,numberOfPaths,0.0);
         if(baseIndices[i][0]!=-1){
             size_t vertex=baseIndices[i][0];
             auto pOutFactor=single_node_cut_factors_[vertex][1]->get_factor();
             const std::unordered_map<size_t,size_t>& baseIDsToIndices=pOutFactor->getBaseIDsToIndices();
             //auto iterT=baseIDsToIndices.find(pInstance->getTerminalNode());
             //assert(iterT!=baseIDsToIndices.end());
             //size_t tIndex=iterT->second;
             lapInput.add_assignment(i,lapInput.no_assignment,pInstance->parameters.getOutputCost());
             pairs[i][numberOfPaths]=pInstance->parameters.getOutputCost();

             //auto iter=pInstance->getMyGraph().forwardNeighborsBegin(vertex);
             size_t j=0;
             for (size_t j = 0; j < numberOfPaths; ++j) {
                 if(baseIndices[j][1]!=-1){

                      size_t vertex2=baseIndices[j][1];
                      assert(vertex2!=pInstance->getTerminalNode());
                     //if(vertex2!=pInstance->getTerminalNode()){
                         auto it=baseIDsToIndices.find(vertex2);
                         //assert(it!=baseIDsToIndices.end());
                         if(it!=baseIDsToIndices.end()){
                             size_t index=it->second;
                             double cost=baseGraph.getForwardEdgeCost(vertex,index);
                             cost+=cummulativeCosts[i][j];
                             lapInput.add_assignment(i,j,cost);
                             pairs[i][j]=cost;
                             if(i==j) origValue+=cost;

                         }
                     //}


                 }
                 else if(j==i){
                     origValue+=pInstance->parameters.getOutputCost();
                 }
             }
         }
         else{
             lapInput.add_assignment(i,lapInput.no_assignment,0.0);
             pairs[i][numberOfPaths]=0;
             origValue+=pInstance->parameters.getInputCost();
         }
         if(baseIndices[i][1]!=-1){
             pairs[numberOfPaths][i]=pInstance->parameters.getInputCost();
             lapInput.add_assignment(lapInput.no_assignment,i,pInstance->parameters.getInputCost());
         }
         else{
             pairs[numberOfPaths][i]=0;
              lapInput.add_assignment(lapInput.no_assignment,i,0.0);

         }
    }
    return origValue;
}


template<class SNC_FACTOR>
void LdpPrimalHeuristics<SNC_FACTOR>::evaluateTimeCuts(){
    const LdpDirectedGraph& liftedGraph=pInstance->getMyGraphLifted();
    const LdpDirectedGraph& baseGraph=pInstance->getMyGraph();



    initCummulativeCosts();

    std::vector<std::array<int,2>> baseIndices=setPathPointers();
    LPMP::linear_assignment_problem_input lapInput;
    double origValue=initLAP(lapInput,baseIndices);


    //TODO: in previous for check if for some i exists a better neighbor than i


    MCF::SSP<long,double> mcf(lapInput.no_mcf_nodes(),lapInput.no_mcf_edges());
    lapInput.initialize_mcf(mcf);

    //  std::cout<<"mcf init"<<std::endl;
    double mcfValue=mcf.solve();

    std::vector<size_t> storeLabeling(numberOfPaths,numberOfPaths+1);
    std::vector<size_t> reverseLabeling(numberOfPaths,numberOfPaths+1);

    std::vector<size_t> permutation(numberOfPaths);
    for (size_t i = 0; i < numberOfPaths; ++i) {
        permutation[i]=i;
    }

    std::vector<size_t> unassignedlabels;
    std::vector<bool> usedLabels(numberOfPaths);
    size_t newNumberOfPaths=numberOfPaths;


   // std::cout<<"orig value "<<origValue<<std::endl;
   // std::cout<<"mcf objective "<<mcfValue<<std::endl;


    for (size_t i = 0; i < mcf.no_edges(); ++i) {
        if(mcf.flow(i)>0.99){
            // int label=mcf.head(i)-numberOfInput;
            assert(lapInput.no_left_nodes==numberOfPaths);
            size_t label=mcf.head(i)-numberOfPaths;
            size_t vertex=mcf.tail(i);
            if(label<numberOfPaths)  usedLabels[label]=1;
            if(vertex<numberOfPaths&&label<numberOfPaths){
                //  std::cout<<"vertex "<<vertex<<", label "<<label<<std::endl;
                storeLabeling[vertex]=label;
                permutation[vertex]=label;
                if(vertex!=label){
                    //std::cout<<"switch"<<std::endl;
                }
                reverseLabeling[label]=vertex;
                if(baseIndices[label][0]==-1){
                  //  std::cout<<"make empty path"<<std::endl;
                    startingVertices[label]=pInstance->getTerminalNode();
                    pointersToPaths[label]=pInstance->getTerminalNode();
                    permutation[label]=numberOfPaths+1;
                }

            }
            else if(label<numberOfPaths&&vertex>=numberOfPaths){
                 //TODO create new path from first path only if this was not a start path
                //unassignedlabels.push_back(label);
                if(baseIndices[label][0]!=-1){   //it will be cut from a path
                    assert(baseIndices[label][1]!=-1);
                    size_t newStartingVertex=baseIndices[label][1];
                    assert(startingVertices.size()==newNumberOfPaths);
                    assert(permutation.size()==newNumberOfPaths);
                    assert(cummulativeCosts.size()==newNumberOfPaths);
                    startingVertices.push_back(newStartingVertex);
                    pointersToPaths.push_back(newStartingVertex);
                    permutation.push_back(label);
                    //neighboringVertices[baseIndices[label][0]]  //TODO later maybe it gets a new neighbor, does not
                    newNumberOfPaths++;



                }

            }
            else if(vertex<numberOfPaths&&label>=numberOfPaths){//this vertex is unassigned, paths needs to be terminated if it was not
                //storeLabeling[vertex]=numberOfPaths+1;
                if(baseIndices[vertex][1]!=-1){   //this path should be terminated
                    //size_t previousDescendant=baseIndices[vertex][1];
                    if(baseIndices[vertex][0]!=-1){
                        permutation[vertex]=numberOfPaths+1;
                    }


                }
                //else permutation stays vertex->vertex

                //TODO create new path from end only if it was not the end of a path
            }
        }
    }


    for (int i = 0; i < numberOfPaths; ++i) {
        size_t label=i;
        if(!usedLabels[i]&&baseIndices[i][1]!=-1&&baseIndices[i][0]!=-1){
            size_t newStartingVertex=baseIndices[label][1];
            //assert(startingVertices.size()==newNumberOfPaths);
            assert(permutation.size()==newNumberOfPaths);
            //assert(cummulativeCosts.size()==newNumberOfPaths);
            startingVertices.push_back(newStartingVertex);
            pointersToPaths.push_back(newStartingVertex);
            permutation.push_back(label);
            //neighboringVertices[baseIndices[label][0]]  //TODO later maybe it gets a new neighbor, does not
            newNumberOfPaths++;


        }

    }

    //Cange cummulative costs according to permutation
    for (size_t i = 0; i < numberOfPaths; ++i) {
        std::vector<double> oldCosts=cummulativeCosts[i];
        std::vector<double> newCosts(newNumberOfPaths);
        assert(permutation.size()==newNumberOfPaths);
        for (int j = 0; j < newNumberOfPaths; ++j) {
            if(permutation[j]<numberOfPaths){
                newCosts[j]=oldCosts[permutation[j]];
            }
            else{
                newCosts[j]=0;
            }

        }
        cummulativeCosts[i]=newCosts;
    }
    for (size_t i = numberOfPaths; i < newNumberOfPaths; ++i) {
        std::vector<double> newCosts(newNumberOfPaths);
        cummulativeCosts.push_back(newCosts);
    }



    //Change list of neighbors and labels according to permutation
    for (size_t i = 0; i < newNumberOfPaths; ++i) {
        if(permutation[i]!=i){

            if(permutation[i]<numberOfPaths){   //connect already existing path
                assert(baseIndices[permutation[i]][1]!=-1);
                size_t secondVertex=baseIndices[permutation[i]][1];
                if(i<numberOfPaths){  //connect to already existing path. Otherwise, starting vertices have already been updated
                    assert(baseIndices[i][0]!=-1);
                    size_t firstVertex=baseIndices[i][0];
                    neighboringVertices[firstVertex]=secondVertex;
                    assert(vertexLabels[firstVertex]==i+1);
                }
                while(secondVertex!=pInstance->getTerminalNode()){
                    vertexLabels[secondVertex]=i+1;
                    secondVertex=neighboringVertices[secondVertex];
                }
            }
            else{ //terminate path if it had a part before current time. Otherwise removin from starting vertices has already been done before.
                assert(i<numberOfPaths);
                if(baseIndices[i][0]!=-1){
                    size_t firstVertex=baseIndices[i][0];
                    neighboringVertices[firstVertex]=pInstance->getTerminalNode();
                }
            }
        }
    }

 numberOfPaths=newNumberOfPaths;
   // currentPrimal();

    //TODO update pointers to paths


}

}
#endif // LDP_PRIMAL_HEURISTICS_HXX
