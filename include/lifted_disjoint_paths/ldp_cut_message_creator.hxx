#ifndef LDP_CUT_MESSAGE_CREATOR_HXX
#define LDP_CUT_MESSAGE_CREATOR_HXX

#include<stdlib.h>
#include<vector>
namespace LPMP {

template <class CUT_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
struct LdpCutMessageInputs{

    void init(CUT_FACTOR* myCutFactor,SINGLE_NODE_CUT_FACTOR_CONT* sncFactor,size_t index);
    std::vector<size_t> _nodeIndicesInCut;
    std::vector<size_t> _nodeIndicesInSnc;
    size_t _nodeIndexOfLiftedEdge;
    bool containsLifted;

};


template <class CUT_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void LdpCutMessageInputs<CUT_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::init(CUT_FACTOR* myCutFactor,SINGLE_NODE_CUT_FACTOR_CONT* snc,size_t i) {
    _nodeIndicesInCut;
    _nodeIndicesInSnc;



    _nodeIndexOfLiftedEdge=0;

    size_t singleVertex=0;
    const auto & outputVertices=myCutFactor->getOutputVertices();
    const auto & inputVertices=myCutFactor->getInputVertices();
    const auto & cutGraph=myCutFactor->getCutGraph();


    containsLifted=false;
    auto * sncFactor=snc->get_factor();
    bool sncIsOut=sncFactor->isNodeOutFlow();


    auto iterCut=cutGraph.forwardNeighborsBegin(0);
    auto iterEnd=cutGraph.forwardNeighborsEnd(0);

    if(!sncIsOut){
        iterCut=cutGraph.backwardNeighborsBegin(i);
        iterEnd=cutGraph.backwardNeighborsEnd(i);
        singleVertex=outputVertices.at(i);
    }
    else{
       iterCut=cutGraph.forwardNeighborsBegin(i);
       iterEnd=cutGraph.forwardNeighborsEnd(i);
        singleVertex=inputVertices.at(i);
    }

    assert(singleVertex==sncFactor->nodeID);

    auto iterSnc=sncFactor->getBaseIDs().begin();
    auto iterSncEnd=sncFactor->getBaseIDs().end();
    if(sncIsOut&&singleVertex==myCutFactor->getLiftedInputVertex()){
        containsLifted=true;
        _nodeIndexOfLiftedEdge=sncFactor->getLiftedIDToOrder(myCutFactor->getLiftedOutputVertex());
        //TODO lifted exists etc
    }
    else if(!sncIsOut&&singleVertex==myCutFactor->getLiftedOutputVertex()){
        containsLifted=true;
        _nodeIndexOfLiftedEdge=sncFactor->getLiftedIDToOrder(myCutFactor->getLiftedInputVertex());
    }

    size_t sncCounter=0;
    size_t cutCounter=0;
    while(iterCut!=iterEnd&&iterSnc!=iterSncEnd){

        size_t vertexInCut=0;
        if(sncIsOut){
            vertexInCut=outputVertices.at(iterCut->head);
        }
        else{
            vertexInCut=inputVertices.at(iterCut->head);
        }
        if(vertexInCut<*iterSnc){
            iterCut++;
            cutCounter++;
        }
        else if(*iterSnc<vertexInCut){
            iterSnc++;
            sncCounter++;
        }
        else {
            assert(*iterSnc==vertexInCut);
            _nodeIndicesInCut.push_back(cutCounter);
            _nodeIndicesInSnc.push_back(sncCounter);
            //snc->updateEdgeCost(-iterCut->second,sncCounter,false);
            iterCut++;
            cutCounter++;
            iterSnc++;
            sncCounter++;
        }
    }

}

}

#endif // LDP_CUT_MESSAGE_CREATOR_HXX
