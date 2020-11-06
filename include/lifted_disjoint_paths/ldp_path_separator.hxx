#ifndef LDP_PATH_SEPARATOR_HXX
#define LDP_PATH_SEPARATOR_HXX

#include"lifted_disjoint_paths/ldp_instance.hxx"
#include"ldp_min_marginals_extractor.hxx"
#include"ldp_path_factor.hxx"

namespace LPMP {

template <class PATH_FACTOR, class PATH_SNC_MESSAGE, class SINGLE_NODE_CUT_FACTOR,class LPFMC>
class ldp_path_separator {

public:
    ldp_path_separator(std::vector<PATH_FACTOR*>& _factorContainer, std::vector<PATH_SNC_MESSAGE*>&  _messageContainer, std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>& _sncFactorContainer, const lifted_disjoint_paths::LdpInstance * _pInstance, ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR>& _mmExtractor,LPFMC *_lp):
    factorContainer(_factorContainer),
    messageContainer(_messageContainer),
    pInstance(_pInstance),
    sncFactorContainer(_sncFactorContainer),
     mmExtractor(_mmExtractor),
      lp_(_lp)
    {

        numberOfVertices=sncFactorContainer.size();
        isInQueue=std::vector<char>(numberOfVertices);
        predInQueue=std::vector<size_t>(numberOfVertices,std::numeric_limits<size_t>::max());
        predInQueueIsLifted=std::vector<char>(numberOfVertices);
    }

    void separatePathInequalities(size_t maxConstraints);
    void createPathFactor(const size_t& lv1,const size_t& lv2,const size_t& bv1,const size_t& bv2,bool isLifted);  //lifted edge vertices and the connecting base edge vertices
    std::list<std::pair<size_t,bool>> findShortestPath(const size_t& firstVertex,const size_t& lastVertex);


private:
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
inline std::list<std::pair<size_t,bool>> ldp_path_separator< PATH_FACTOR,  PATH_SNC_MESSAGE, SINGLE_NODE_CUT_FACTOR,LPFMC>::findShortestPath(const size_t& firstVertex,const size_t& lastVertex){
   std::list<std::pair<size_t,bool>> shortestPath;
   //BSF should work best
   //I need a map of used edges
   std::list<size_t> queue; //vertex and is lifted (between vertex and its predecessor in queue)

   const auto& vg=pInstance->getVertexGroups();
   size_t firstIndex=vg.getGroupIndex(firstVertex);
   size_t lastIndex=vg.getGroupIndex(lastVertex);
   queue.push_back(firstVertex) ;  //is lifted for the first does not matter
   bool pathFound=false;
   std::vector<size_t> toDeleteFromQueue;
   toDeleteFromQueue.push_back(firstIndex);

   //TODO I need pointers to predecessors
   while(!queue.empty()&&!pathFound){
       size_t& vertex=queue.front();

       std::vector<std::pair<size_t,bool>>& neighbors=usedEdges[vertex];

       for(size_t i=0;i<neighbors.size();i++){
           std::pair<size_t,bool> p=neighbors[i];
           size_t neighborOfVertex=p.first;
           bool isLifted=p.second;
           if(neighborOfVertex==lastVertex){
               predInQueue[lastVertex]=vertex;
               predInQueueIsLifted[lastVertex]=isLifted;
               isInQueue[lastVertex]=1; //To be cleared

               pathFound=true;
               break;
           }

           size_t nodeTimeIndex=vg.getGroupIndex(neighborOfVertex);
           if(nodeTimeIndex<lastIndex&&!isInQueue[neighborOfVertex]){
               queue.push_back(neighborOfVertex);
               isInQueue[neighborOfVertex]=1;
               predInQueue[neighborOfVertex]=vertex;
               predInQueueIsLifted[neighborOfVertex]=isLifted;
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
    std::list<std::pair<size_t,bool>> beginning=findShortestPath(lv1,bv1);
    std::list<std::pair<size_t,bool>> ending=findShortestPath(bv2,lv2);
    std::vector<size_t> pathVertices(beginning.size()+ending.size());
    std::vector<char> liftedEdgesIndices(beginning.size()+ending.size());
    auto iter=beginning.begin();
    auto end=beginning.end();
    size_t counter=0;
    for (;iter!=end;iter++) {
        size_t vertex=iter->first;
        bool isLifted=iter->second;
        pathVertices[counter]=vertex;
        liftedEdgesIndices[counter]=isLifted;
        counter++;
    }
    pathVertices[counter]=bv1;
    liftedEdgesIndices[counter]=isLifted;
    counter++;
    //Careful with the bridge edge!
    iter=ending.begin();
    end=ending.end();
    for(;iter!=end;iter++){
        size_t vertex=iter->first;
        bool isLifted=iter->second;
        pathVertices[counter]=vertex;
        liftedEdgesIndices[counter]=isLifted;
        counter++;
    }
    pathVertices[counter]=lv2;
    liftedEdgesIndices[counter]=true; //this relates to the big lifted edge connecting the first and the last path vertex

    std::vector<double> costs(pathVertices.size(),0); //last member is the lifted edge cost
    auto* pathFactor = lp_->template add_factor<PATH_FACTOR>(pathVertices,costs,liftedEdgesIndices);
    factorContainer.push_back(pathFactor);






}

template <class PATH_FACTOR, class PATH_SNC_MESSAGE, class SINGLE_NODE_CUT_FACTOR,class LPFMC>
inline void ldp_path_separator< PATH_FACTOR,  PATH_SNC_MESSAGE, SINGLE_NODE_CUT_FACTOR,LPFMC>::separatePathInequalities(size_t maxConstraints){
    mmExtractor.initMinMarginals();
    baseMM=mmExtractor.getBaseEdgesMinMarginals();
    liftedMM=mmExtractor.getLiftedEdgesMinMarginals();



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

    for(size_t i=0;i<edgesToSort.size();i++){
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
            while(iterDescPred!=endDescPred&&iterDescV2!=endDescV2){
                if(*iterDescV2<*iterDescPred){ //exists desc of v2 not contained in desc of pred
                    const size_t& descV2=*iterDescV2;
                    auto f=positiveLifted[pred].find(descV2);
                    if(f!=positiveLifted[pred].end()){  //Contradicting lifted edge exists!
                        createPathFactor(pred,descV2,vertex1,vertex2); //TODO first just put to a queue (list of vertices and information if the edges are lifted) and then select the best
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









}



}



#endif // LDP_PATH_SEPARATOR_HXX
