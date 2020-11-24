#ifndef LDP_CUT_FACTOR_SEPARATOR_HXX
#define LDP_CUT_FACTOR_SEPARATOR_HXX

#include"lifted_disjoint_paths/ldp_instance.hxx"
#include"ldp_min_marginals_extractor.hxx"
#include<list>

namespace LPMP {


template <class CUT_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
class LdpCutSeparator{


public:
    LdpCutSeparator(const lifted_disjoint_paths::LdpInstance * _pInstance, ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR_CONT>& _mmExtractor):

    pInstance(_pInstance),
     mmExtractor(_mmExtractor)

    {
        numberOfVertices=pInstance->getNumberOfVertices()-2;

    }

    void separateCutInequalities(size_t maxConstraints);

    std::priority_queue<std::pair<double,CUT_FACTOR*>>& getPriorityQueue(){
        return pQueue;
    }


    //LdpPathMessageInputs getMessageInputsToPathFactor(PATH_FACTOR* myPathFactor,SINGLE_NODE_CUT_FACTOR_CONT* sncFactor,size_t index)const ;
    void clearPriorityQueue();


private:
    void connectEdge(const size_t& v,const size_t& w);
    void createCut(size_t v1,size_t v2,double cost);

    const lifted_disjoint_paths::LdpInstance * pInstance;
    ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR_CONT>& mmExtractor;
    size_t numberOfVertices;


    std::vector<std::map<size_t,double>> baseEdgesWithCosts;
    std::vector<std::map<size_t,double>> liftedEdgesWithCosts;
    std::vector<std::list<size_t>> predecessors;
    std::vector<std::list<size_t>> descendants;
    std::vector<std::vector<char>> isConnected;
    std::priority_queue<std::pair<double,CUT_FACTOR*>> pQueue;
     std::vector<std::list<size_t>> candidateLifted;

};

template <class CUT_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void LdpCutSeparator<CUT_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::clearPriorityQueue() {
    while(!pQueue.empty()){
        std::pair<double,CUT_FACTOR*> p=pQueue.top();
        delete p.second;
        p.second=nullptr;
        pQueue.pop();
    }

}

template  <class CUT_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void LdpCutSeparator<CUT_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::connectEdge(const size_t& v,const size_t& w){
    assert(v<numberOfVertices);
    assert(w<numberOfVertices);
    const LdpDirectedGraph & baseGraph=pInstance->getMyGraph();
    for(auto& pred: predecessors[v]){
        assert(pred<numberOfVertices);
        auto  origDescendants=descendants[pred].begin();
        auto  newDescendants=descendants[w].begin();
        auto  endOrig=descendants[pred].end();
        auto  endNew=descendants[w].end();
        auto  itBase=baseGraph.forwardNeighborsBegin(pred);
        size_t baseCounter=0;
        auto  baseEnd=baseGraph.forwardNeighborsEnd(pred);
        //TODO update is connected



        while(newDescendants!=endNew){
            while(itBase!=baseEnd&&itBase->first<*newDescendants){
                itBase++;
                baseCounter++;
            }
            if(origDescendants==endOrig||*origDescendants>*newDescendants){
                descendants[pred].insert(origDescendants,(*newDescendants));
                predecessors[*newDescendants].push_back(pred);
                if(itBase!=baseEnd&&itBase->first==*newDescendants){
                    assert(baseCounter<isConnected[pred].size());
                    isConnected[pred][baseCounter]=1;
                    baseCounter++;
                    itBase++;
                }
                newDescendants++;
            }
            else if(*newDescendants>*origDescendants){
                origDescendants++;
            }
            else{
                assert(*origDescendants=*newDescendants);
                origDescendants++;
                newDescendants++;
            }
        }
    }
}


template  <class CUT_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void LdpCutSeparator<CUT_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::createCut(size_t v1,size_t v2,double cost){

    double lCost=liftedEdgesWithCosts[v1][v2];  //TODO maybe obtain without map?
    std::map<size_t,std::map<size_t,double>> cutEdges;
    const LdpDirectedGraph & baseGraph=pInstance->getMyGraph();

    double improvementValue=std::min(abs(lCost),cost);

    auto descV1Iter=descendants[v1].begin();
    auto descV1end=descendants[v1].end();
    for (;descV1Iter!=descV1end;descV1Iter++) {
        auto descV1SecondIter=descV1Iter;
        size_t d=*descV1Iter;
        const auto* it=baseGraph.forwardNeighborsBegin(d);
        const auto* end=baseGraph.forwardNeighborsEnd(d);
        while(it!=end){
            if(descV1SecondIter==descV1end||it->first<*descV1SecondIter){
                size_t d2=it->first;
                if(pInstance->isReachable(d2,v2)){
                    assert(baseEdgesWithCosts[d][d2]>=cost-eps);
                    cutEdges[d][d2]=0;
                    assert(d<numberOfVertices);
                    assert(d2<numberOfVertices);
                }
                it++;
            }
            else if(it->first>*descV1SecondIter){
                descV1SecondIter++;
            }
            else {
                assert(*descV1SecondIter==it->first);
                it++;
                descV1SecondIter++;
            }

        }


    }

//    for(auto& d:descendants[v1]){              //or break if d is in used vertices for cuts and w is reachable from d
//        const auto* it=baseGraph.forwardNeighborsBegin(d);
//        const auto* end=baseGraph.forwardNeighborsEnd(d);
//        for(;it!=end;it++){               //can be probably a linear iteration
//            size_t d2=it->first;

//            //TODO fix this! How to effectively check if d2 is descendant of v1?
//            if(!isConnected[v1][d2]&&pInstance->isReachable(d2,v2)){  //candidate cut edge (d,d2), leaves component and w is reachable
////                if(baseEdgesWithCosts[d][d2]<cost-eps){
////                    std::cout<<"wrong cut edge "<<d<<", "<<d2<<", value: "<<cost<<"is connected "<<isConnected[d][d2]<<std::endl;
////                }
//                assert(baseEdgesWithCosts[d][d2]>=cost-eps);
//                cutEdges[d][d2]=0;
//                assert(d<numberOfVertices);
//                assert(d2<numberOfVertices);
//            }
//        }

//    }
    // ldp_cut_factor* pCutF=new ldp_cut_factor(v1,v2,lCost,cutEdges); //TODO use this after solving upddate cost in snc
    CUT_FACTOR* pCutF=new CUT_FACTOR(v1,v2,0.0,cutEdges);
    pQueue.push(std::pair(improvementValue,pCutF));
}



template  <class CUT_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void LdpCutSeparator<CUT_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::separateCutInequalities(size_t maxConstraints){
    std::cout<<"separate cuts "<<std::endl;
    mmExtractor.initMinMarginals();
    baseEdgesWithCosts=mmExtractor.getBaseEdgesMinMarginals();
    liftedEdgesWithCosts=mmExtractor.getLiftedEdgesMinMarginals();



    std::vector<std::tuple<double,size_t,size_t>> edgesToSort;
    descendants= std::vector<std::list<size_t>> (numberOfVertices);
    predecessors= std::vector<std::list<size_t>> (numberOfVertices);


    const LdpDirectedGraph & baseGraph=pInstance->getMyGraph();


    isConnected=std::vector<std::vector<char>>(numberOfVertices);
    for (size_t i = 0; i < numberOfVertices; ++i) {
        size_t s=baseGraph.getNumberOfEdgesFromVertex(i);
        isConnected[i]=std::vector<char>(s,0);
    }

    //Structure for connecting negative (and small positive?) edges
    for(size_t node=0;node<predecessors.size();node++){
        predecessors[node].push_back(node);
        descendants[node].push_back(node);
    }

    std::cout<<"desc and pred init done"<<std::endl;


    //list of base edges to be sorted
    for(size_t i=0;i<baseEdgesWithCosts.size();i++){
    //for(auto it=baseEdgesWithCosts.begin();it!=baseEdgesWithCosts.end();it++){
        size_t vertex=i;
        std::map<size_t,double>& neighbors=baseEdgesWithCosts[i];
        size_t neighborsCounter=0;
        for(auto it2=neighbors.begin();it2!=neighbors.end();it2++, neighborsCounter++){
            size_t w=it2->first;
            double cost=it2->second;
            assert(baseGraph.getForwardEdgeVertex(vertex,neighborsCounter)==w);
            //std::tuple<double,size_t,size_t> t(cost,v,w)
            if(vertex!=pInstance->getSourceNode()&&w!=pInstance->getTerminalNode()){
                if(cost<eps){
                    connectEdge(vertex,w);
                }
                else{

                    edgesToSort.push_back(std::tuple<double,size_t,size_t>(cost,vertex,neighborsCounter));
                }
            }
        }
    }

    std::sort(edgesToSort.begin(),edgesToSort.end());

    std::cout<<"edges to sort sorted"<<std::endl;


    //Select candidate lifted edges: negative and disconnected

    candidateLifted=std::vector<std::list<size_t>> (numberOfVertices);  //pair size_t,double instead of size_t?
    size_t nrClosedNodes=0;
    std::vector<char> closedNodes(numberOfVertices);
    for (size_t i=0;i<numberOfVertices;i++) {
        const std::map<size_t,double>& neighbors=liftedEdgesWithCosts.at(i);
        //std::map<size_t,double> neighborsToKeep;
        auto itLifted=neighbors.begin();
        auto itDesc=descendants[i].begin();
        while(itLifted!=neighbors.end()){
            if(itDesc==descendants[i].end()||*itDesc>itLifted->first){
                if(itLifted->second<eps) candidateLifted[i].push_back(itLifted->first);
                itLifted++;
            }
            else if(*itDesc<itLifted->first){
                itDesc++;
            }
            else{
                assert(*itDesc==itLifted->first);
                itDesc++;
                itLifted++;
            }
        }
        if(candidateLifted[i].empty()){
            nrClosedNodes++;
        }

    }

    std::cout<<"candidate lifted obtained "<<std::endl;




    std::cout<<"number of vertices "<<numberOfVertices<<std::endl;
    size_t i=0;
    while(nrClosedNodes<numberOfVertices&&i<edgesToSort.size()){
        size_t v=std::get<1>(edgesToSort[i]);
        size_t index=std::get<2>(edgesToSort[i]);
        size_t w=baseGraph.getForwardEdgeVertex(v,index);
       // double cost=std::get<0>(edgesToSort[i]);

        assert(v<isConnected.size());
        assert(index<isConnected[v].size());
        if(isConnected[v][index]){
            i++;
            continue;
        }
        double cost=std::get<0>(edgesToSort[i]);
        assert(cost>0);

        std::cout<<v<<", "<<w<<":"<<cost<<std::endl;
        std::cout<<"number of closed "<<nrClosedNodes<<std::endl;

      //  std::tuple<double,size_t,size_t> bestLiftedEdge;  //lb improvement, cost, vertices
      //  double bestLiftedCost=0;
      //  std::cout<<"compute"<<cost<<std::endl;

        for(auto& pred: predecessors[v]){
            std::cout<<"pred "<<pred<<std::endl;
            if(candidateLifted[pred].empty()) continue;
            std::list<size_t> edgesToKeep;
            auto  liftedIt=candidateLifted[pred].begin();
            auto  liftedEnd=candidateLifted[pred].end();
            auto  newDescIt=descendants[w].begin();
            auto  newDescEnd=descendants[w].end();
            auto  oldDescIt=descendants[pred].begin();
            auto  oldDescEnd=descendants[pred].end();

            while(newDescIt!=newDescEnd&&liftedIt!=liftedEnd){
                std::cout<<"new desc, lifted "<<(*newDescIt)<<", "<<(*liftedIt)<<std::endl;
                while(liftedIt!=liftedEnd&&*liftedIt<*newDescIt) liftedIt++;
                if(oldDescIt==oldDescEnd||*oldDescIt>*newDescIt){
                    std::cout<<"new desc "<<std::endl;
                    if(liftedIt!=liftedEnd&&*liftedIt==*newDescIt){
                        createCut(pred,*liftedIt,cost);
                        liftedIt=candidateLifted[pred].erase(liftedIt);
                    }
                    newDescIt++;
                }
                else if(*oldDescIt<*newDescIt){
                    oldDescIt++;
                }
                else{
                    assert(*oldDescIt==*newDescIt);
                   // createCut(desc,*liftedIt);
                    oldDescIt++;
                    newDescIt++;
                }
            }
            if(candidateLifted[pred].empty()){
                nrClosedNodes++;
            }
        }

        connectEdge(v,w);
        i++;
    }

    std::cout<<"queue with cuts filled"<<std::endl;


}







}

#endif // LDP_CUT_FACTOR_SEPARATOR_HXX
