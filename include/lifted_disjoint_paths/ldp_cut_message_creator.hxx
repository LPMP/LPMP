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

    //bool _sncIsOut;

    _nodeIndexOfLiftedEdge;

    size_t inputVertex=myCutFactor->getInputVertices().at(i);
    const auto & outputVertices=myCutFactor->getOutputVertices();
    const auto & cutGraph=myCutFactor->getCutGraph();


    containsLifted=false;
    auto * sncFactorOut=snc->get_factor();
    assert(inputVertex==sncFactorOut->nodeID);
    auto iterCut=cutGraph.forwardNeighborsBegin(i);
    auto iterEnd=cutGraph.forwardNeighborsEnd(i);

    auto iterSnc=sncFactorOut->getBaseIDs().begin();
    auto iterSncEnd=sncFactorOut->getBaseIDs().end();
    if(inputVertex==myCutFactor->getLiftedInputVertex()){
        containsLifted=true;
        _nodeIndexOfLiftedEdge=sncFactorOut->getLiftedIDToOrder(myCutFactor->getLiftedOutputVertex());
        //TODO lifted exists etc
    }
    size_t sncCounter=0;
    size_t cutCounter=0;
    while(iterCut!=iterEnd&&iterSnc!=iterSncEnd){
        size_t vertexInCut=outputVertices.at(iterCut->head);
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
