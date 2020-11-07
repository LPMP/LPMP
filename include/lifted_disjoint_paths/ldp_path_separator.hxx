#ifndef LDP_PATH_SEPARATOR_HXX
#define LDP_PATH_SEPARATOR_HXX

#include"lifted_disjoint_paths/ldp_instance.hxx"
#include"ldp_min_marginals_extractor.hxx"
#include"ldp_path_factor.hxx"

namespace LPMP {


struct LdpPathMessageInputs{
    std::vector<size_t> edgeIndexInPath;  //Mostly one, two for the first and the last path vertices
    std::vector<size_t> vertexIndexInSnc;
    std::vector<char> isLifted;

};


template <class PATH_FACTOR, class PATH_SNC_MESSAGE, class SINGLE_NODE_CUT_FACTOR,class LPFMC>
class ldp_path_separator {

public:
    ldp_path_separator(std::vector<PATH_FACTOR*>& _factorContainer, std::vector<PATH_SNC_MESSAGE*>&  _messageContainer, std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>& _sncFactorContainer, const lifted_disjoint_paths::LdpInstance * _pInstance, ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR>& _mmExtractor,LPFMC *_lp):
    factorContainer(_factorContainer),
    messageContainer(_messageContainer),
    pInstance(_pInstance),
    sncFactorContainer(_sncFactorContainer),
     mmExtractor(_mmExtractor)

    {
        lp_=_lp;

        numberOfVertices=sncFactorContainer.size();
        isInQueue=std::vector<char>(numberOfVertices);
        predInQueue=std::vector<size_t>(numberOfVertices,std::numeric_limits<size_t>::max());
        predInQueueIsLifted=std::vector<char>(numberOfVertices);
    }

    size_t separatePathInequalities(size_t maxConstraints);

    LdpPathMessageInputs getMessageInputsToPathFactor(PATH_FACTOR* pathFactor,SINGLE_NODE_CUT_FACTOR* sncFactor,size_t index,bool isOut)const ;


private:

    void createPathFactor(const size_t& lv1,const size_t& lv2,const size_t& bv1,const size_t& bv2,bool isLifted);  //lifted edge vertices and the connecting base edge vertices
    std::list<std::pair<size_t,bool>> findShortestPath(const size_t& firstVertex,const size_t& lastVertex);


std::vector<PATH_FACTOR*>& factorContainer;
std::vector<PATH_SNC_MESSAGE*>&  messageContainer;
const lifted_disjoint_paths::LdpInstance * pInstance;
std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>& sncFactorContainer;
ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR>& mmExtractor;
size_t numberOfVertices;
LPFMC *lp_;

std::vector<std::map<size_t,double>> baseMM;
std::vector<std::map<size_t,double>> liftedMM;
std::vector<std::list<size_t>> predecessors;  //Can I use list? Maybe yes, just predecessors will not be sorted!
std::vector<std::list<size_t>> descendants;
 std::vector<std::vector<std::pair<size_t,bool>>> usedEdges;
std::vector<char> isInQueue;
std::vector<size_t> predInQueue;
std::vector<char> predInQueueIsLifted;

};

template <class PATH_FACTOR, class PATH_SNC_MESSAGE, class SINGLE_NODE_CUT_FACTOR,class LPFMC>
inline LdpPathMessageInputs ldp_path_separator< PATH_FACTOR,  PATH_SNC_MESSAGE, SINGLE_NODE_CUT_FACTOR,LPFMC>::getMessageInputsToPathFactor(PATH_FACTOR* pathFactor,SINGLE_NODE_CUT_FACTOR* sncFactor,size_t index,bool isOut)const {
    auto * myPathFactor=pathFactor->get_factor();

    size_t numberOfEdges=myPathFactor->getNumberOfEdges();
    assert(index<numberOfEdges);

    const std::vector<size_t>& pathVertices=myPathFactor->getListOfVertices();
    const std::vector<char>& liftedInfo=myPathFactor->getLiftedInfo();

    std::vector<size_t> indicesInSnc;
    std::vector<char> isLiftedForMessage;
    std::vector<size_t> edgeIndicesInPath;
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
        indicesInSnc[1]=pFactorInFactor->getLiftedIDToOrder(pathVertices.back());

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

    LdpPathMessageInputs messageInputs={edgeIndicesInPath,indicesInSnc,isLiftedForMessage};
    return messageInputs;

}


template <class PATH_FACTOR, class PATH_SNC_MESSAGE, class SINGLE_NODE_CUT_FACTOR,class LPFMC>
inline std::list<std::pair<size_t,bool>> ldp_path_separator< PATH_FACTOR,  PATH_SNC_MESSAGE, SINGLE_NODE_CUT_FACTOR,LPFMC>::findShortestPath(const size_t& firstVertex,const size_t& lastVertex){
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



template <class PATH_FACTOR, class PATH_SNC_MESSAGE, class SINGLE_NODE_CUT_FACTOR,class LPFMC>
inline void ldp_path_separator< PATH_FACTOR,  PATH_SNC_MESSAGE, SINGLE_NODE_CUT_FACTOR,LPFMC>::createPathFactor(const size_t& lv1, const size_t& lv2, const size_t &bv1, const size_t &bv2,bool isLifted){

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

   // std::cout<<lv2<<", "<<std::endl;
    //std::cout<<"lifted "<<isLifted<<std::endl;

    std::vector<double> costs(pathVertices.size(),0); //last member is the lifted edge cost
    assert(pathVertices[0]==lv1);
    assert(pathVertices.back()==lv2);

    auto* pathFactor = lp_->template add_factor<PATH_FACTOR>(pathVertices,costs,liftedEdgesIndices);
    factorContainer.push_back(pathFactor);

    assert(lv1<numberOfVertices);
    auto * pFactorOut=sncFactorContainer[lv1][1];
    auto * pFactorOutFactor=pFactorOut->get_factor();

    std::vector<size_t> indicesInSnc={0,0};
    std::vector<char> isLiftedForFirstMessage={0,1};
    std::vector<size_t> edgeIndicesInPath={0,numberOfVerticesInPath-1};

    //size_t secondPathVertexIndex;
    if(liftedEdgesIndices[0]){
        isLiftedForFirstMessage[0]=1;
        indicesInSnc[0]=pFactorOutFactor->getLiftedIDToOrder(pathVertices[1]);
    }
    else{
        isLiftedForFirstMessage[0]=0;
        indicesInSnc[0]=pFactorOutFactor->getBaseIDToOrder(pathVertices[1]);
    }
    indicesInSnc[1]=pFactorOutFactor->getLiftedIDToOrder(lv2);

   // std::cout<<"edge indices size "<<edgeIndicesInPath.size()<<std::endl;
    ldp_snc_path_message myMessage(edgeIndicesInPath,indicesInSnc,isLiftedForFirstMessage);
    //auto * message=lp_->template add_message<PATH_SNC_MESSAGE>(pathFactor,pFactorOut,edgeIndicesInPath,indicesInSnc,isLiftedForFirstMessage);
    auto * message=lp_->template add_message<PATH_SNC_MESSAGE>(pathFactor,pFactorOut,myMessage);
  //  std::cout<<"message start added, dimension "<<message->get_message()->dimension<<std::endl;
    messageContainer.push_back(message);
  //  std::cout<<"message in container, indices dimension "<<edgeIndicesInPath.size()<<std::endl;

    assert(lv2<numberOfVertices);

    auto * pFactorIn=sncFactorContainer[lv2][0];
    auto * pFactorInFactor=pFactorIn->get_factor();

   // std::cout<<"pointers set "<<std::endl;
    indicesInSnc={0,0};
    isLiftedForFirstMessage={0,1};
    edgeIndicesInPath={numberOfVerticesInPath-2,numberOfVerticesInPath-1};

   // std::cout<<"arrays set "<<std::endl;

    if(liftedEdgesIndices[numberOfVerticesInPath-2]){
      //  std::cout<<"last edge lifted"<<std::endl;
        isLiftedForFirstMessage[0]=1;
        indicesInSnc[0]=pFactorInFactor->getLiftedIDToOrder(pathVertices[numberOfVerticesInPath-2]);
    }
    else{
      //  std::cout<<"last edge base"<<std::endl;
        isLiftedForFirstMessage[0]=0;
        indicesInSnc[0]=pFactorInFactor->getBaseIDToOrder(pathVertices[numberOfVerticesInPath-2]);
    }
   // std::cout<<"big edge is lifted "<<std::endl;
    indicesInSnc[1]=pFactorInFactor->getLiftedIDToOrder(lv1);

    //std::cout<<"try to add second message" <<std::endl;
    //auto * message2=lp_->template add_message<PATH_SNC_MESSAGE>(pathFactor,pFactorIn,edgeIndicesInPath,indicesInSnc,isLiftedForFirstMessage);
    ldp_snc_path_message myMessage2(edgeIndicesInPath,indicesInSnc,isLiftedForFirstMessage);
    auto * message2=lp_->template add_message<PATH_SNC_MESSAGE>(pathFactor,pFactorIn,myMessage2);
    messageContainer.push_back(message2);
   // std::cout<<"message end added "<<std::endl;

    indicesInSnc=std::vector<size_t>(1);
    isLiftedForFirstMessage=std::vector<char>(1);
    edgeIndicesInPath=std::vector<size_t>(1);


    for(size_t i=1;i<numberOfVerticesInPath-1;i++){
        pFactorIn=sncFactorContainer[pathVertices[i]][0];
        pFactorInFactor=pFactorIn->get_factor();
        edgeIndicesInPath[0]=i;
        if(liftedEdgesIndices[i-1]){
            isLiftedForFirstMessage[0]=1;
            indicesInSnc[0]=pFactorInFactor->getLiftedIDToOrder(pathVertices[i-1]);
        }
        else{
            isLiftedForFirstMessage[0]=0;
            indicesInSnc[0]=pFactorInFactor->getBaseIDToOrder(pathVertices[i-1]);
        }

        ldp_snc_path_message myMessage3(edgeIndicesInPath,indicesInSnc,isLiftedForFirstMessage);
        auto * message3=lp_->template add_message<PATH_SNC_MESSAGE>(pathFactor,pFactorIn,myMessage3);
        messageContainer.push_back(message3);

//        message=lp_->template add_message<PATH_SNC_MESSAGE>(pathFactor,pFactorIn,edgeIndicesInPath,indicesInSnc,isLiftedForFirstMessage);
//        messageContainer.push_back(message);

     //   std::cout<<"message added "<<std::endl;
        pFactorOut=sncFactorContainer[pathVertices[i]][1];
        pFactorOutFactor=pFactorOut->get_factor();

        if(liftedEdgesIndices[i]){
            isLiftedForFirstMessage[0]=1;
            indicesInSnc[0]=pFactorOutFactor->getLiftedIDToOrder(pathVertices[i+1]);
        }
        else{
            isLiftedForFirstMessage[0]=0;
            indicesInSnc[0]=pFactorOutFactor->getBaseIDToOrder(pathVertices[i+1]);
        }

        ldp_snc_path_message myMessage4(edgeIndicesInPath,indicesInSnc,isLiftedForFirstMessage);
        auto * message4=lp_->template add_message<PATH_SNC_MESSAGE>(pathFactor,pFactorIn,myMessage4);
        messageContainer.push_back(message4);

//        message2=lp_->template add_message<PATH_SNC_MESSAGE>(pathFactor,pFactorOut,edgeIndicesInPath,indicesInSnc,isLiftedForFirstMessage);
//        messageContainer.push_back(message2);

    }




    std::cout<<"message container size "<<messageContainer.size()<<std::endl;





}

template <class PATH_FACTOR, class PATH_SNC_MESSAGE, class SINGLE_NODE_CUT_FACTOR,class LPFMC>
inline size_t ldp_path_separator< PATH_FACTOR,  PATH_SNC_MESSAGE, SINGLE_NODE_CUT_FACTOR,LPFMC>::separatePathInequalities(size_t maxConstraints){
    mmExtractor.initMinMarginals();
    baseMM=mmExtractor.getBaseEdgesMinMarginals();
    liftedMM=mmExtractor.getLiftedEdgesMinMarginals();


    size_t constraintsCounter=0;

    std::vector<std::tuple<double,size_t,size_t,bool>> edgesToSort; //contains negative base and lifted edges: cost,vertex1,vertex2,isLifted

    for(size_t i=0;i<numberOfVertices;i++){
        auto iter=baseMM[i].begin();
        auto end=baseMM[i].end();
        for(;iter!=end;iter++){
            if(iter->second<0){
                edgesToSort.push_back(std::tuple(iter->second,i,iter->first,false));
            }
        }
    }

    std::vector<std::map<size_t,double>> positiveLifted(numberOfVertices);
    for(size_t i=0;i<numberOfVertices;i++){
        auto iter=liftedMM[i].begin();
        auto end=liftedMM[i].end();
        for(;iter!=end;iter++){
            if(iter->second<0){
                edgesToSort.push_back(std::tuple(iter->second,i,iter->first,true));
            }
            else{
                positiveLifted[i][iter->first]=iter->second;
            }
        }
    }

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

    for(size_t i=0;i<edgesToSort.size()&&constraintsCounter<maxConstraints;i++){
        std::tuple<double,size_t,size_t,bool>& edge=edgesToSort[i];
        size_t& vertex1=std::get<1>(edge);
        size_t& vertex2=std::get<2>(edge);
        bool isLifted=std::get<3>(edge);

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
                        createPathFactor(pred,descV2,vertex1,vertex2,isLifted); //TODO first just put to a queue (list of vertices and information if the edges are lifted) and then select the best
                        constraintsCounter++;
                    }

                    descendants[pred].insert(iterDescPred,descV2);
                    predecessors[descV2].push_back(pred);
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

    std::cout<<"path constraints "<<constraintsCounter<<std::endl;
    return  constraintsCounter;










}



}



#endif // LDP_PATH_SEPARATOR_HXX
