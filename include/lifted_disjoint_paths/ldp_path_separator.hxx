#ifndef LDP_PATH_SEPARATOR_HXX
#define LDP_PATH_SEPARATOR_HXX

#include"lifted_disjoint_paths/ldp_instance.hxx"
#include"ldp_min_marginals_extractor.hxx"
#include"ldp_path_factor.hxx"

namespace LPMP {

template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
struct LdpPathMessageInputs{

    void init(PATH_FACTOR* myPathFactor,SINGLE_NODE_CUT_FACTOR_CONT* sncFactor,size_t index);

    std::vector<size_t> edgeIndicesInPath;  //Mostly one, two for the first and the last path vertices
    std::vector<size_t> indicesInSnc;
    std::vector<char> isLiftedForMessage;

};


template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void LdpPathMessageInputs<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::init(PATH_FACTOR* myPathFactor,SINGLE_NODE_CUT_FACTOR_CONT* sncFactor,size_t index) {
  //  auto * myPathFactor=pathFactor->get_factor();

    bool isOut=sncFactor->get_factor()->isNodeOutFlow();

    size_t numberOfEdges=myPathFactor->getNumberOfEdges();
    assert(index<numberOfEdges);

    const std::vector<size_t>& pathVertices=myPathFactor->getListOfVertices();
    assert(pathVertices.at(index)==sncFactor->get_factor()->nodeID);
    const std::vector<char>& liftedInfo=myPathFactor->getLiftedInfo();


    if(index==0){
        assert(isOut);
        auto * pFactorOutFactor=sncFactor->get_factor();
        indicesInSnc={0,0};
        isLiftedForMessage={liftedInfo.at(0),1};
        edgeIndicesInPath={0,numberOfEdges-1};


        //size_t secondPathVertexIndex;
        if(liftedInfo.at(0)){
            indicesInSnc[0]=pFactorOutFactor->getLiftedIDToOrder(pathVertices.at(1));
        }
        else{
            indicesInSnc[0]=pFactorOutFactor->getBaseIDToOrder(pathVertices.at(1));
        }
        indicesInSnc[1]=pFactorOutFactor->getLiftedIDToOrder(pathVertices.back());

    }
    else if(index==numberOfEdges-1){
        assert(!isOut);
        auto * pFactorInFactor=sncFactor->get_factor();
        indicesInSnc={0,0};
        isLiftedForMessage={liftedInfo.at(numberOfEdges-2),1};
        edgeIndicesInPath={numberOfEdges-2,numberOfEdges-1};

        if(liftedInfo.at(numberOfEdges-2)){
            indicesInSnc[0]=pFactorInFactor->getLiftedIDToOrder(pathVertices.at(numberOfEdges-2));
        }
        else{
            indicesInSnc[0]=pFactorInFactor->getBaseIDToOrder(pathVertices.at(numberOfEdges-2));
        }
        indicesInSnc[1]=pFactorInFactor->getLiftedIDToOrder(pathVertices[0]);

    }
    else if(isOut){
        auto * pFactorOutFactor=sncFactor->get_factor();
        indicesInSnc={0};
        isLiftedForMessage={liftedInfo.at(index)};
        edgeIndicesInPath={index};
        if(liftedInfo.at(index)){
            indicesInSnc[0]=pFactorOutFactor->getLiftedIDToOrder(pathVertices.at(index+1));
        }
        else{
            indicesInSnc[0]=pFactorOutFactor->getBaseIDToOrder(pathVertices.at(index+1));
        }


    }
    else {
        auto * pFactorInFactor=sncFactor->get_factor();
        indicesInSnc={0};
        isLiftedForMessage={liftedInfo.at(index-1)};
        edgeIndicesInPath={index-1};
        if(liftedInfo.at(index-1)){
            indicesInSnc[0]=pFactorInFactor->getLiftedIDToOrder(pathVertices.at(index-1));
        }
        else{
            indicesInSnc[0]=pFactorInFactor->getBaseIDToOrder(pathVertices.at(index-1));
        }
    }

//    LdpPathMessageInputs messageInputs={edgeIndicesInPath,indicesInSnc,isLiftedForMessage};
//    return messageInputs;

}


template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT> //PATH_FACTOR is ldp_path_factor, SINGLE_NODE_CUT_FACTOR_CONT is the container wrapper
class ldp_path_separator {

public:
    ldp_path_separator(const lifted_disjoint_paths::LdpInstance * _pInstance, ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR_CONT>& _mmExtractor):

    pInstance(_pInstance),
     mmExtractor(_mmExtractor)

    {
        numberOfVertices=pInstance->getNumberOfVertices()-2;
        isInQueue=std::vector<char>(numberOfVertices);
        predInQueue=std::vector<size_t>(numberOfVertices,std::numeric_limits<size_t>::max());
        predInQueueIsLifted=std::vector<char>(numberOfVertices);
        maxTimeGap=std::max(pInstance->getGapLifted(),pInstance->getGapBase());
    }

    void separatePathInequalities(size_t maxConstraints);

    std::priority_queue<std::pair<double,PATH_FACTOR*>>& getPriorityQueue(){
        return pQueue;
    }


    //LdpPathMessageInputs getMessageInputsToPathFactor(PATH_FACTOR* myPathFactor,SINGLE_NODE_CUT_FACTOR_CONT* sncFactor,size_t index)const ;
    void clearPriorityQueue();


private:

    PATH_FACTOR* createPathFactor(const size_t& lv1,const size_t& lv2,const size_t& bv1,const size_t& bv2,bool isLifted);  //lifted edge vertices and the connecting base edge vertices
    std::list<std::pair<size_t,bool>> findShortestPath(const size_t& firstVertex,const size_t& lastVertex);

const lifted_disjoint_paths::LdpInstance * pInstance;
ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR_CONT>& mmExtractor;
size_t numberOfVertices;


std::vector<std::map<size_t,double>> baseMM;
std::vector<std::map<size_t,double>> liftedMM;
std::vector<std::list<size_t>> predecessors;  //Can I use list? Maybe yes, just predecessors will not be sorted!
std::vector<std::list<size_t>> descendants;
 std::vector<std::vector<std::pair<size_t,bool>>> usedEdges;
std::vector<char> isInQueue;
std::vector<size_t> predInQueue;
std::vector<char> predInQueueIsLifted;
std::priority_queue<std::pair<double,PATH_FACTOR*>> pQueue;
 size_t maxTimeGap;

};

template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void ldp_path_separator<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::clearPriorityQueue() {
    while(!pQueue.empty()){
        std::pair<double,PATH_FACTOR*> p=pQueue.top();
        delete p.second;
        p.second=nullptr;
        pQueue.pop();
    }

}
  //  auto * myPathFactor=pathFactor->get_factor();



template <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline std::list<std::pair<size_t,bool>> ldp_path_separator<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::findShortestPath(const size_t& firstVertex,const size_t& lastVertex){
   assert(firstVertex!=lastVertex);

   std::list<std::pair<size_t,bool>> shortestPath;
   //BSF should work best
   //I need a map of used edges
   std::list<size_t> queue; //vertex and is lifted (between vertex and its predecessor in queue)

  // std::cout<<"finding shortest path between "<<firstVertex<<", "<<lastVertex<<std::endl;
   const auto& vg=pInstance->getVertexGroups();
   size_t firstIndex=vg.getGroupIndex(firstVertex);
   size_t lastIndex=vg.getGroupIndex(lastVertex);
 //  std::cout<<"time of last vertex "<<lastIndex<<std::endl;
   queue.push_back(firstVertex) ;  //is lifted for the first does not matter
   bool pathFound=false;
   std::vector<size_t> toDeleteFromQueue;
   toDeleteFromQueue.push_back(firstIndex);

   //TODO I need pointers to predecessors
   while(!queue.empty()&&!pathFound){
       size_t& vertex=queue.front();

       std::vector<std::pair<size_t,bool>>& neighbors=usedEdges[vertex];

      // std::cout<<"vertex in queue "<<vertex<<", neighbors "<<std::endl;
       for(size_t i=0;i<neighbors.size();i++){
           std::pair<size_t,bool> p=neighbors[i];
           size_t neighborOfVertex=p.first;
           bool isLifted=p.second;
          // std::cout<<neighborOfVertex<<", "<<std::endl;
           if(neighborOfVertex==lastVertex){
               //std::cout<<"is last vertex"<<std::endl;
               predInQueue[lastVertex]=vertex;
               predInQueueIsLifted[lastVertex]=isLifted;
               isInQueue[lastVertex]=1; //To be cleared
               toDeleteFromQueue.push_back(lastVertex);

               pathFound=true;
               break;
           }

           size_t nodeTimeIndex=vg.getGroupIndex(neighborOfVertex);
           if(nodeTimeIndex<lastIndex&&!isInQueue[neighborOfVertex]){
             //  std::cout<<"gets to queue "<<std::endl;
               queue.push_back(neighborOfVertex);
               isInQueue[neighborOfVertex]=1;
               predInQueue[neighborOfVertex]=vertex;
               predInQueueIsLifted[neighborOfVertex]=isLifted;
               toDeleteFromQueue.push_back(neighborOfVertex);
           }
       }
       queue.pop_front();
   }
   assert(pathFound);

   size_t currentVertex=lastVertex;
   while(currentVertex!=firstVertex){  //path contains vertex and info about edge starting in it
       size_t newVertex=predInQueue[currentVertex];
       bool isEdgeLifted=predInQueueIsLifted[currentVertex];
       shortestPath.push_front({newVertex,isEdgeLifted});
       currentVertex=newVertex;
   }

   size_t maxValue=std::numeric_limits<size_t>::max();
   for(auto& v:toDeleteFromQueue){
       isInQueue[v]=0;
       predInQueue[v]=maxValue;
       predInQueueIsLifted[v]=0;
   }
   return shortestPath;
}



template  <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline PATH_FACTOR* ldp_path_separator<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::createPathFactor(const size_t& lv1, const size_t& lv2, const size_t &bv1, const size_t &bv2,bool isLifted){

    std::list<std::pair<size_t,bool>> beginning;
    if(lv1!=bv1) beginning=findShortestPath(lv1,bv1);
    std::list<std::pair<size_t,bool>> ending;
    if(bv2!=lv2) ending=findShortestPath(bv2,lv2);
    if(lv1==bv1&&lv2==bv2) std::cout<<"base covering lifted for path "<<std::endl;
    std::vector<size_t> pathVertices(beginning.size()+ending.size()+2);
    std::vector<char> liftedEdgesIndices(beginning.size()+ending.size()+2);


    size_t numberOfVerticesInPath=pathVertices.size();

    auto iter=beginning.begin();
    auto end=beginning.end();
    size_t counter=0;
   // std::cout<<"path vertices"<<std::endl;
    for (;iter!=end;iter++) {
        size_t vertex=iter->first;
        bool isLiftedEdge=iter->second;
       // std::cout<<vertex<<", "<<std::endl;
      //  std::cout<<"lifted "<<isLiftedEdge<<std::endl;
        pathVertices[counter]=vertex;
        liftedEdgesIndices[counter]=isLiftedEdge;
        counter++;
    }
    pathVertices[counter]=bv1;
    liftedEdgesIndices[counter]=isLifted;
  //  std::cout<<"central edge vertex "<<bv1<<", "<<std::endl;
  //  std::cout<<"lifted "<<isLifted<<std::endl;
    counter++;
    //Careful with the bridge edge!
    if(bv2!=lv2) assert(ending.front().first==bv2);
    iter=ending.begin();
    end=ending.end();
    for(;iter!=end;iter++){
        size_t vertex=iter->first;
        bool isLiftedEdge=iter->second;
      //  std::cout<<vertex<<", "<<std::endl;
      //  std::cout<<"lifted "<<isLiftedEdge<<std::endl;
        pathVertices[counter]=vertex;
        liftedEdgesIndices[counter]=isLiftedEdge;
        counter++;
    }
    pathVertices[counter]=lv2;
    liftedEdgesIndices[counter]=true; //this relates to the big lifted edge connecting the first and the last path vertex

    assert(counter==numberOfVerticesInPath-1);
        assert(pathVertices.front()<pInstance->getNumberOfVertices()-2&&pathVertices.back()<pInstance->getNumberOfVertices());

   // std::cout<<lv2<<", "<<std::endl;
    //std::cout<<"lifted "<<isLifted<<std::endl;

    std::vector<double> costs(pathVertices.size(),0); //last member is the lifted edge cost
    assert(pathVertices[0]==lv1);
    assert(pathVertices.back()==lv2);

    PATH_FACTOR* pPathFactor=new PATH_FACTOR(pathVertices,costs,liftedEdgesIndices,pInstance);
    return pPathFactor;



}

template  <class PATH_FACTOR,class SINGLE_NODE_CUT_FACTOR_CONT>
inline void ldp_path_separator<PATH_FACTOR,SINGLE_NODE_CUT_FACTOR_CONT>::separatePathInequalities(size_t maxConstraints){
    mmExtractor.initMinMarginals();
    baseMM=mmExtractor.getBaseEdgesMinMarginals();
    liftedMM=mmExtractor.getLiftedEdgesMinMarginals();


    size_t constraintsCounter=0;

    std::vector<std::tuple<double,size_t,size_t,bool>> edgesToSort; //contains negative base and lifted edges: cost,vertex1,vertex2,isLifted

    for(size_t i=0;i<numberOfVertices;i++){
        auto iter=baseMM[i].begin();
        auto end=baseMM[i].end();
        for(;iter!=end;iter++){
            if(iter->first<numberOfVertices&&iter->second<-eps){
                edgesToSort.push_back(std::tuple(iter->second,i,iter->first,false));
            }
        }
    }

    std::vector<std::map<size_t,double>> positiveLifted(numberOfVertices);
    for(size_t i=0;i<numberOfVertices;i++){
        auto iter=liftedMM[i].begin();
        auto end=liftedMM[i].end();
        for(;iter!=end;iter++){
            if(iter->second<-eps){
                edgesToSort.push_back(std::tuple(iter->second,i,iter->first,true));
            }
            else if(iter->second>eps){
                positiveLifted[i][iter->first]=iter->second;
            }
        }
    }  maxTimeGap=std::max(pInstance->getGapLifted(),pInstance->getGapBase());

    std::sort(edgesToSort.begin(),edgesToSort.end());

//    std::vector<std::set<size_t>> predecessors(numberOfVertices);  //Can I use list? Maybe yes, just predecessors will not be sorted!
//    std::vector<std::set<size_t>> descendants(numberOfVertices);

//    for(size_t i=0;i<numberOfVertices;i++){
//        predecessors[i].insert(i);
//        descendants[i].insert(i);

//    }

    predecessors=std::vector<std::list<size_t>> (numberOfVertices);  //Can I use list? Maybe yes, just predecessors will not be sorted!
    descendants=std::vector<std::list<size_t>> (numberOfVertices);

    for(size_t i=0;i<numberOfVertices;i++){
        predecessors[i].push_back(i);
        descendants[i].push_back(i);

    }

    usedEdges= std::vector<std::vector<std::pair<size_t,bool>>> (numberOfVertices);


    for(size_t i=0;i<edgesToSort.size();i++){
        std::tuple<double,size_t,size_t,bool>& edge=edgesToSort[i];
        size_t& vertex1=std::get<1>(edge);
        size_t& vertex2=std::get<2>(edge);
        bool isLifted=std::get<3>(edge);
        double edgeCost=std::get<0>(edge);

        auto iterPredV1=predecessors[vertex1].begin();
        auto endPredV1=predecessors[vertex1].end();
        auto iterDescV2=descendants[vertex2].begin();
        auto endDescV2=descendants[vertex2].end();

        for (;iterPredV1!=endPredV1;iterPredV1++) {
            const size_t& pred=*iterPredV1;
            auto iterDescPred=descendants[pred].begin();  //Put descendants of V2 into descendants of pred
            auto endDescPred=descendants[pred].end();
            while(iterDescV2!=endDescV2){
                if(iterDescPred==endDescPred||*iterDescV2<*iterDescPred){ //exists desc of v2 not contained in desc of pred
                    const size_t& descV2=*iterDescV2;
                   // std::cout<<"exists new descendant "<<std::endl;
                    auto f=positiveLifted[pred].find(descV2);
                    if(f!=positiveLifted[pred].end()){  //Contradicting lifted edge exists!
                        if(pred!=vertex1||descV2!=vertex2){
                            PATH_FACTOR* pPathFactor= createPathFactor(pred,descV2,vertex1,vertex2,isLifted); //TODO first just put to a queue (list of vertices and information if the edges are lifted) and then select the best
                            double improvementValue=std::min(abs(edgeCost),f->second);
                            pQueue.push(std::pair(improvementValue,pPathFactor));
                            constraintsCounter++;
                        }
                    }

                    size_t l0=pInstance->getGroupIndex(pred);
                    size_t l1=pInstance->getGroupIndex(descV2);
                    if(l1-l0<=maxTimeGap){
                        descendants[pred].insert(iterDescPred,descV2);
                        predecessors[descV2].push_back(pred);
                    }
                    iterDescV2++;

                }
                else if (*iterDescV2>*iterDescPred) { //not interesting
                    iterDescPred++;
                }
                else{  //not interesting
                    iterDescPred++;
                    iterDescV2++;
                }
            }
        }
        usedEdges[vertex1].push_back(std::pair(vertex2,isLifted));

        //TODO add predecessors and descendanta here
    }

    mmExtractor.clearMinMarginals();

    std::cout<<"candidate constraints "<<constraintsCounter<<std::endl;







}



}



#endif // LDP_PATH_SEPARATOR_HXX
