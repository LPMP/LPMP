#ifndef LDP_GAEC_HXX
#define LDP_GAEC_HXX

#include<stdlib.h>
#include<map>
#include<vector>
#include<lifted_disjoint_paths/ldp_directed_graph.hxx>
#include "ldp_instance.hxx"


namespace LPMP {


template <class SNC_FACTOR>
class LdpGreedyAdditive{

public:
    LdpGreedyAdditive(std::vector<std::array<SNC_FACTOR*,2>>* _pSNCFactors, const lifted_disjoint_paths::LdpInstance * _pInstance):
        vg(pInstance->getVertexGroups())
    {
        pSNCFactors=_pSNCFactors;
        pInstance=_pInstance;

        numberOfVertices=pSNCFactors->size();
        t=numberOfVertices+2;
        s=numberOfVertices+1;

        costs=std::vector<std::vector<double>>(numberOfVertices);
        //liftedCosts=std::vector<std::map<size_t,double>>(numberOfVertices);
        liftedCostsFull=std::vector<ShiftedVector<double>>(numberOfVertices);
        timeStamps=std::vector<std::vector<size_t>>(numberOfVertices);

        startingVertices=std::vector<char>(numberOfVertices,1);
      //  endingVertices=std::vector<char>(numberOfVertices,1);

        neighbors=std::vector<size_t>(numberOfVertices,t);

        inOutCost=pInstance->parameters.getInputCost()+pInstance->parameters.getOutputCost();


        const std::vector<std::array<SNC_FACTOR*,2>>& single_node_cut_factors_=*pSNCFactors;

        for (size_t i = 0; i < numberOfVertices; ++i) { //init costs, lifted costs and time stamps, costs in forward direction only

            auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
            costs[i]=pOutFactor->getBaseCosts();
            timeStamps[i]=std::vector<size_t>(costs[i].size());

            std::vector<double> lc=pOutFactor->getLiftedCosts();
            std::vector<size_t> liftedIDs=pOutFactor->getLiftedCosts();

            liftedCostsFull[i]=ShiftedVector<double>(getMinForwardNeighbor(i),getMaxForwardNeighbor(i),0);
            assert(lc.size()==liftedIDs.size());
            for (size_t j = 0; j < liftedIDs.size(); ++j) {
                size_t vertex2=liftedIDs[j];
                //liftedCosts[i][vertex2]=lc[j];
                liftedCostsFull[i][vertex2]=lc[j];
            }

        }

        for (size_t i = 0; i < numberOfVertices; ++i) { //update costs with opposite direction

            auto pInFactor=single_node_cut_factors_[i][0]->get_factor();
            const std::vector<double>& baseCosts=pInFactor->getBaseCosts();
            const std::vector<size_t>& baseIDs=pInFactor->getBaseIDs();

            for (int j = 0; j < baseCosts.size(); ++j) {
                size_t vertexID=baseIDs[j];
                if(vertexID>=numberOfVertices) continue;
                auto pOutFactor=single_node_cut_factors_[vertexID][1]->get_factor();

                size_t index=pOutFactor->getBaseIDToOrder(i);
                assert(index<costs[vertexID].size());
                costs[vertexID][index]+=baseCosts[j];

            }

            const std::vector<double>& localLiftedCosts=pInFactor->getLiftedCosts();
            const std::vector<size_t>& liftedIDs=pInFactor->getLiftedIDs();

            for (int j = 0; j < liftedIDs.size(); ++j) {
                size_t vertexID=liftedIDs[j];
                //liftedCosts[vertexID][i]+=localLiftedCosts[j];
                liftedCostsFull[vertexID][i]+=localLiftedCosts[j];
            }
        }

//        for (int i = 0; i < costs.size(); ++i) {
//            auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
//            const std::vector<size_t>& liftedIDs=pOutFactor->getLiftedIDs();
//            const std::vector<size_t>& baseIDs=pOutFactor->getBaseIDs();

//            for(size_t j=0;j<baseIDs.size();j++){  //add lifted costs to base costs
//                    costs[i][j]+=liftedCostsFull[i][baseIDs[j]];
//            }
//        }

        for (size_t i = 0; i < costs.size(); ++i) {
            auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
            const std::vector<size_t>& baseIDs=pOutFactor->getBaseIDs();

            assert(baseIDs.size()==costs[i].size());
            for (size_t j = 0; j < costs[i].size(); ++j) {
                size_t vertex2=baseIDs[j];
                if(vertex2>=numberOfVertices) continue;
                size_t timeGap=getTimeGap(i,vertex2);
                double cost=costs[i][j]+liftedCostsFull[i][vertex2];
                if(cost<inOutCost){
                    Edge e={i,j,cost,timeGap,0};
                    Q.push(e);
                }
            }
        }



    }


    void runGAEC(){

        while(!Q.empty()){
            Edge e=Q.top();
            size_t secondVertex=getVertexToIndex(e.firstVertex);
            size_t timeStamp=timeStamps[e.firstVertex][e.neighborIndex];
            if(startingVertices[secondVertex]&&neighbors[e.firstVertex]==t&&timeStamp==e.timeStamp){
                startingVertices[secondVertex]=false;
                //endingVertices[e.firstVertex]=false;
                neighbors[e.firstVertex]=secondVertex;
                updateCostOutgoing(e.firstVertex,secondVertex);
                updateCostIncomming(e.firstVertex,secondVertex);
            }

            Q.pop();
        }
    }


private:
    const std::vector<size_t>& getStartingVertices()const{
        return startingFinal;
    }

    const std::vector<size_t>& getNeighbors()const{
        return neighbors;
    }

    void updateCostOutgoing(size_t i, size_t j){
        assert(j>i);
        ShiftedVector<double>& neighborsOfFirst=liftedCostsFull[i];
        ShiftedVector<double>& neighborsOfSecond=liftedCostsFull[j];
        for (size_t l = neighborsOfSecond.getMinVertex(); l <= neighborsOfFirst.getMaxVertex(); ++l) {
            neighborsOfSecond[l]+=neighborsOfFirst[l];
        }
        //neighborsOfFirst=ShiftedVector<double>();

        auto pOutFactor=(*pSNCFactors)[j][1]->get_factor();
        const std::vector<size_t>& baseIDs=pOutFactor->getBaseIDs();

        assert(baseIDs.size()==costs[j].size());
        for (size_t k = 0; k < costs[j].size(); ++k) {
            size_t vertex2=baseIDs[k];
            if(vertex2==t) continue;
            size_t timeGap=getTimeGap(j,vertex2);
            double cost=costs[j][k]+neighborsOfSecond[vertex2];
            timeStamps[j][k]++;
            if(cost<inOutCost){
                Edge e={j,k,cost,timeGap,timeStamps[j][k]};
                Q.push(e);
            }
        }
    }


    void updateCostIncomming(size_t i, size_t j){
        assert(j>i);

        size_t minCommonVertex=getMinBackwardNeighbor(j);
        size_t maxCommonVertex=getMaxBackwardNeighbor(i);

        if(maxCommonVertex!=s){
            for (size_t l = minCommonVertex; l <= maxCommonVertex; ++l) {
                liftedCostsFull[l][i]+=liftedCostsFull[l][j];
            }
        }

        auto pInFactor=(*pSNCFactors)[i][0]->get_factor();
        const std::vector<size_t>& baseIDs=pInFactor->getBaseIDs();

        //assert(baseIDs.size()==costs[j].size());
        for (size_t k = 0; k < baseIDs.size(); ++k) {
            size_t vertex2=baseIDs[k];
            if(vertex2==s) continue;
            auto pOutFactor=(*pSNCFactors)[vertex2][1]->get_factor();
            size_t index=pOutFactor->getBaseIDToOrder(i);
            size_t timeGap=getTimeGap(vertex2,i);
            double cost=costs[vertex2][index]+liftedCostsFull[vertex2][i];
            timeStamps[vertex2][index]++;
            if(cost<inOutCost){
                Edge e={vertex2,index,cost,timeGap,timeStamps[vertex2][index]};
                Q.push(e);
            }
        }
    }

    size_t getMinForwardNeighbor(size_t vertex){
        size_t time=vg.getGroupIndex(vertex);
        if(time==vg.getMaxTime()){
            return t;
        }
        else{
            return vg.getMinVertexInTime(time+1);
        }
    }

    size_t getMaxForwardNeighbor(size_t vertex){
        size_t time=vg.getGroupIndex(vertex);
        size_t shiftedTime=time+pInstance->parameters.getMaxTimeGapComplete();
        if(time==vg.getMaxTime()){
            return t;
        }
        else if(shiftedTime>=vg.getMaxTime()){
            return vg.getMaxVertex();
        }
        else{
            return vg.getMaxVertexInTime(shiftedTime);
        }
    }

    size_t getMaxBackwardNeighbor(size_t vertex){
        size_t time=vg.getGroupIndex(vertex);
        if(time==1){
            return s;
        }
        else{
            return vg.getMinVertexInTime(time-1);
        }
    }

    size_t getMinBackwardNeighbor(size_t vertex){
        size_t time=vg.getGroupIndex(vertex);
        //size_t shiftedTime=time+pInstance->parameters.getMaxTimeGapComplete();
        if(time==1){
            return s;
        }
        else if(time<=pInstance->parameters.getMaxTimeGapComplete()){
            return 0;
        }
        else{
            return vg.getMaxVertexInTime(time-pInstance->parameters.getMaxTimeGapComplete());
        }
    }



    size_t getTimeGap(size_t vertex1,size_t vertex2){
        size_t l0=pInstance->getGroupIndex(vertex1);
        size_t l1=pInstance->getGroupIndex(vertex2);

        assert(l1>l0);
        return l1-l0;
    }

    size_t getTimeForIndex(size_t vertex1,size_t index){
        assert(vertex1<numberOfVertices);
        auto pOutFactor=(*pSNCFactors)[vertex1][1]->get_factor();
        const std::vector<size_t>& baseIDs=pOutFactor->getBaseIDs();
        assert(index<baseIDs.size());
        return getTimeGap(vertex1,baseIDs[index]);
    }

    size_t getVertexToIndex(size_t vertex1,size_t index){
        assert(vertex1<numberOfVertices);
        auto pOutFactor=(*pSNCFactors)[vertex1][1]->get_factor();
        const std::vector<size_t>& baseIDs=pOutFactor->getBaseIDs();
        assert(index<baseIDs.size());
        return baseIDs.at(index);
    }

    void finalizeResults(){

        for (std::size_t graph_node = 0; graph_node < numberOfVertices; ++graph_node) {
            (*pSNCFactors)[graph_node][0]->get_factor()->setNoBaseEdgeActive();
            (*pSNCFactors)[graph_node][1]->get_factor()->setNoBaseEdgeActive();
        }

        for (size_t i = 0; i < numberOfVertices; ++i) {
            if(startingVertices[i]){
               // if(neighbors[i]!=t){
                    startingFinal.push_back(startingVertices[i]);
                    auto pInFactor=(*pSNCFactors)[i][0]->get_factor();
                    pInFactor->setBaseEdgeActiveWithID(s);
               // }
                //else{//is isolated

                //}

            }
        }

        for (int i = 0; i < numberOfVertices; ++i) {
            size_t neighbor=neighbors[i];
            auto pOutFactor=(*pSNCFactors)[i][1]->get_factor();
            pOutFactor->setBaseEdgeActiveWithID(neighbor);
            if(neighbor!=t){
                auto pInFactor=(*pSNCFactors)[neighbor][0]->get_factor();
                pInFactor->setBaseEdgeActiveWithID(i);
            }

        }
    }

    struct Edge
    {
         size_t firstVertex;
         size_t neighborIndex;
         double cost;
         size_t timeGap;
         size_t timeStamp;

         bool operator <(Edge const& other) const
         {
             return cost/double(timeGap) < other.cost/double(other.timeGap);
         }

    };

    size_t numberOfVertices;
    size_t t;
    size_t s;
    //std::vector<std::map<size_t, double>> costs;   //TODO replace with LdpDirectedGraph
   // std::vector<std::map<size_t, double>> liftedCosts; //TODO replace with LdpDirectedGraph

    std::vector<ShiftedVector<double>> liftedCostsFull;

    std::vector<std::vector<double>> costs;   //TODO replace with LdpDirectedGraph
    //std::vector<std::vector<double>> liftedCosts; //TODO replace with LdpDirectedGraph

    std::vector<std::vector<size_t>> timeStamps;

    std::vector<char> startingVertices; //char for representing bool: is starting, to check if can be connected in forward direction
    std::vector<size_t> startingFinal;
   // std::vector<char> endingVertices;  //to check if it can be connected in a backward direction

    std::vector<size_t> neighbors;

    std::vector<std::array<SNC_FACTOR*,2>>* pSNCFactors;

    const lifted_disjoint_paths::LdpInstance* pInstance;

     std::priority_queue<Edge> Q;

     const VertexGroups<>& vg;

     double inOutCost;

     //std::vector<std::array<size_t,2>> boundariesOfTimeLayers;  //max vertex and min vertex in each time layer


};



//LdpGreedyAdditive(std::vector<std::array<SNC_FACTOR*,2>>* _pSNCFactors, const lifted_disjoint_paths::LdpInstance * _pInstance){
//    pSNCFactors=_pSNCFactors;
//    pInstance=_pInstance;

//    numberOfVertices=pSNCFactors->size();
//    t=numberOfVertices+2;

//    costs=std::vector<std::vector<double>>(numberOfVertices);
//    liftedCosts=std::vector<std::vector<double>>(numberOfVertices);
//    timeStamps=std::vector<std::vector<size_t>>(numberOfVertices);

//    startingVertices=std::vector<char>(numberOfVertices,1);
//    endingVertices=std::vector<char>(numberOfVertices,1);

//    neighbors=std::vector<size_t>(numberOfVertices,t);

//    const std::vector<std::array<SNC_FACTOR*,2>>& single_node_cut_factors_=*pSNCFactors;

//    for (size_t i = 0; i < numberOfVertices; ++i) { //init costs, lifted costs and time stamps, costs in forward direction only

//        auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
//        costs[i]=pOutFactor->getBaseCosts();
//        timeStamps[i]=std::vector<size_t>(costs[i].size());
//        liftedCosts[i]=pOutFactor->getLiftedCosts();
//    }

//    for (size_t i = 0; i < numberOfVertices; ++i) { //update costs with opposite direction

//        auto pInFactor=single_node_cut_factors_[i][0]->get_factor();
//        const std::vector<double>& baseCosts=pInFactor->getBaseCosts();
//        const std::vector<size_t>& baseIDs=pInFactor->getBaseIDs();

//        for (int j = 0; j < baseCosts.size(); ++j) {
//            size_t vertexID=baseIDs[j];
//            auto pOutFactor=single_node_cut_factors_[vertexID][1]->get_factor();

//            size_t index=pOutFactor->getBaseIDToOrder(i);
//            assert(index<costs[vertexID].size());
//            costs[vertexID][index]+=baseCosts[j];

//        }

//        const std::vector<double>& localLiftedCosts=pInFactor->getLiftedCosts();
//        const std::vector<size_t>& liftedIDs=pInFactor->getLiftedIDs();

//        for (int j = 0; j < liftedIDs.size(); ++j) {
//            size_t vertexID=liftedIDs[j];
//            auto pOutFactor=single_node_cut_factors_[vertexID][1]->get_factor();

//            size_t indexLifted=pOutFactor->getLiftedIDToOrder(i);
//            assert(indexLifted<liftedCosts[vertexID].size());
//            liftedCosts[vertexID][indexLifted]+=localLiftedCosts[j];

//        }
//    }

//    for (int i = 0; i < costs.size(); ++i) {
//        auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
//        const std::vector<size_t>& liftedIDs=pOutFactor->getLiftedIDs();
//        const std::vector<size_t>& baseIDs=pOutFactor->getBaseIDs();


//        size_t baseCounter=0;
//        size_t liftedCounter=0;

//        while(liftedCounter<liftedIDs.size()&&baseCounter<baseIDs.size()){  //add lifted costs to base costs
//            while(liftedCounter<liftedIDs.size()&&liftedIDs[liftedCounter]<baseIDs[baseCounter]) liftedCounter++;
//            if(liftedIDs[liftedCounter]==baseIDs[baseCounter]){
//                costs[i][baseCounter]+=liftedCosts[i][liftedCounter];
//            }
//            baseCounter++;
//        }
//    }

//    for (size_t i = 0; i < costs.size(); ++i) {
//        auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
//        const std::vector<size_t>& baseIDs=pOutFactor->getBaseIDs();

//        assert(baseIDs.size()==costs[i].size());
//        for (size_t j = 0; j < costs[i].size(); ++j) {
//            size_t vertex2=baseIDs[j];
//            size_t timeGap=getTimeGap(i,vertex2);
//            Edge e={i,j,costs[i][j],timeGap,0};
//            Q.push(e);
//        }
//    }



//}



//LdpGreedyAdditive(std::vector<std::array<SNC_FACTOR*,2>>* _pSNCFactors, const lifted_disjoint_paths::LdpInstance * _pInstance):
//    vg(pInstance->getVertexGroups())
//{
//    pSNCFactors=_pSNCFactors;
//    pInstance=_pInstance;

//    numberOfVertices=pSNCFactors->size();
//    t=numberOfVertices+2;
//    s=numberOfVertices+1;

//    costs=std::vector<std::vector<double>>(numberOfVertices);
//    liftedCosts=std::vector<std::map<size_t,double>>(numberOfVertices);
//    timeStamps=std::vector<std::vector<size_t>>(numberOfVertices);

//    startingVertices=std::vector<char>(numberOfVertices,1);
//    endingVertices=std::vector<char>(numberOfVertices,1);

//    neighbors=std::vector<size_t>(numberOfVertices,t);


//    const std::vector<std::array<SNC_FACTOR*,2>>& single_node_cut_factors_=*pSNCFactors;

//    for (size_t i = 0; i < numberOfVertices; ++i) { //init costs, lifted costs and time stamps, costs in forward direction only

//        auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
//        costs[i]=pOutFactor->getBaseCosts();
//        timeStamps[i]=std::vector<size_t>(costs[i].size());

//        std::vector<double> lc=pOutFactor->getLiftedCosts();
//        std::vector<size_t> liftedIDs=pOutFactor->getLiftedCosts();

//        assert(lc.size()==liftedIDs.size());
//        for (size_t j = 0; j < liftedIDs.size(); ++j) {
//            size_t vertex2=liftedIDs[j];
//            liftedCosts[i][vertex2]=lc[j];
//        }

//    }

//    for (size_t i = 0; i < numberOfVertices; ++i) { //update costs with opposite direction

//        auto pInFactor=single_node_cut_factors_[i][0]->get_factor();
//        const std::vector<double>& baseCosts=pInFactor->getBaseCosts();
//        const std::vector<size_t>& baseIDs=pInFactor->getBaseIDs();

//        for (int j = 0; j < baseCosts.size(); ++j) {
//            size_t vertexID=baseIDs[j];
//            auto pOutFactor=single_node_cut_factors_[vertexID][1]->get_factor();

//            size_t index=pOutFactor->getBaseIDToOrder(i);
//            assert(index<costs[vertexID].size());
//            costs[vertexID][index]+=baseCosts[j];

//        }

//        const std::vector<double>& localLiftedCosts=pInFactor->getLiftedCosts();
//        const std::vector<size_t>& liftedIDs=pInFactor->getLiftedIDs();

//        for (int j = 0; j < liftedIDs.size(); ++j) {
//            size_t vertexID=liftedIDs[j];
//            liftedCosts[vertexID][i]+=localLiftedCosts[j];
//        }
//    }

//    for (int i = 0; i < costs.size(); ++i) {
//        auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
//        const std::vector<size_t>& liftedIDs=pOutFactor->getLiftedIDs();
//        const std::vector<size_t>& baseIDs=pOutFactor->getBaseIDs();


//        size_t baseCounter=0;
//        auto pLifted=liftedCosts[i].begin();

//        while(pLifted!=liftedCosts[i].end()&&baseCounter<baseIDs.size()){  //add lifted costs to base costs
//            while(pLifted!=liftedCosts[i].end()&&pLifted->first<baseIDs[baseCounter]) pLifted++;
//            if(pLifted->first==baseIDs[baseCounter]){
//                costs[i][baseCounter]+=pLifted->second;
//            }
//            baseCounter++;
//        }
//    }

//    for (size_t i = 0; i < costs.size(); ++i) {
//        auto pOutFactor=single_node_cut_factors_[i][1]->get_factor();
//        const std::vector<size_t>& baseIDs=pOutFactor->getBaseIDs();

//        assert(baseIDs.size()==costs[i].size());
//        for (size_t j = 0; j < costs[i].size(); ++j) {
//            size_t vertex2=baseIDs[j];
//            size_t timeGap=getTimeGap(i,vertex2);
//            Edge e={i,j,costs[i][j],timeGap,0};
//            Q.push(e);
//        }
//    }



//}



}
#endif // LDP_GAEC_HXX
