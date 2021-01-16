#include"lifted_disjoint_paths/ldp_instance.hxx"


namespace LPMP{
namespace lifted_disjoint_paths {






LdpInstance::LdpInstance(LdpParameters<> &configParameters, CompleteStructure<>& cs):
    parameters(configParameters)
{




    pCompleteGraph=&cs.myCompleteGraph;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    myGraphLifted=cs.myCompleteGraph;
    vertexGroups=cs.getVertexGroups();

    if(parameters.isMustCutMissing()){
        initCanJoinStructure(cs.myCompleteGraph);
    }

    s_=myGraphLifted.getNumberOfVertices();
    t_=s_+1;
    myGraph=LdpDirectedGraph(cs.myCompleteGraph,parameters.getInputCost(),parameters.getOutputCost());
    vertexScore=cs.verticesScore;



    if(parameters.isUseAdaptiveThreshold()){
        if(diagnostics()) std::cout<<"using adaptive"<<std::endl;
       //initAdaptiveThresholds(&cs.completeScore,nullptr);
       initAdaptiveThresholds(&myGraph,&myGraphLifted);
    }
    else{
        baseThreshold=parameters.getBaseUpperThreshold();
        positiveLiftedThreshold=parameters.getPositiveThresholdLifted();
        negativeLiftedThreshold=parameters.getNegativeThresholdLifted();
    }

    numberOfVertices=myGraph.getNumberOfVertices();

    init();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const std::chrono::steady_clock::time_point& csBegin=cs.getContructorBegin();
    parameters.getControlOutput()<< "Time of instance constructor = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << " seconds" << std::endl;
    parameters.getControlOutput()<< "Time of instance constructor and complete structure = " << std::chrono::duration_cast<std::chrono::seconds> (end - csBegin).count() << " seconds" << std::endl;
    parameters.writeControlOutput();

    parameters.getControlOutput()<< "Edge file name " << cs.getEdgeFileName()<< std::endl;
    parameters.writeControlOutput();




}

LdpInstance::LdpInstance(LdpParameters<>& configParameters,LdpBatchProcess& BP):
    parameters(configParameters){

    pCompleteGraph=&BP.getMyCompleteGraph();
    parameters.getControlOutput()<< "Vertices in the complete graph "<<pCompleteGraph->getNumberOfVertices()<< std::endl;
    parameters.writeControlOutput();

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    vertexScore=BP.getVerticesScore();

    BP.createLocalVG(vertexGroups);


    myGraphLifted=BP.getMyCompleteGraph();
    s_=myGraphLifted.getNumberOfVertices();
    t_=s_+1;
    myGraph=LdpDirectedGraph(BP.getMyCompleteGraph(),parameters.getInputCost(),parameters.getOutputCost());

    if(parameters.isMustCutMissing()){
        initCanJoinStructure(myGraphLifted);
    }

    if(parameters.isUseAdaptiveThreshold()){
        if(diagnostics()) std::cout<<"using adaptive"<<std::endl;

       initAdaptiveThresholds(&myGraph,&myGraphLifted);
    }
    else{
        baseThreshold=parameters.getBaseUpperThreshold();
        positiveLiftedThreshold=parameters.getPositiveThresholdLifted();
        negativeLiftedThreshold=parameters.getNegativeThresholdLifted();
    }
    numberOfVertices=myGraph.getNumberOfVertices();

    init();


    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

     parameters.getControlOutput() << "Time of instance constructor = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << " seconds" << std::endl;
    const std::chrono::steady_clock::time_point& bpBegin=BP.getContructorBegin();
    parameters.getControlOutput() << "Time of instance constructor and batch process = " << std::chrono::duration_cast<std::chrono::seconds> (end - bpBegin).count() << " seconds" << std::endl;
    parameters.writeControlOutput();


}









void LdpInstance::initAdaptiveThresholds(const LdpDirectedGraph* pBaseGraph,const LdpDirectedGraph* pLiftedGraph){

    const LdpDirectedGraph *pToSort;

    if(pLiftedGraph!=nullptr){
        pToSort=pLiftedGraph;

    }
    else{
        pToSort=pBaseGraph;
    }
    assert(parameters.getPositiveThresholdLifted()>=0&&parameters.getPositiveThresholdLifted()<1);
    assert(parameters.getNegativeThresholdLifted()>=0&&parameters.getNegativeThresholdLifted()<1);
    assert(parameters.getBaseUpperThreshold()>=0&&parameters.getBaseUpperThreshold()<1);



    std::set<double> lowestCosts;

    size_t numberOfNegative=0;
    for (int i = 0; i < pToSort->getNumberOfVertices(); ++i) {
        auto iter=pToSort->forwardNeighborsBegin(i);
        for (;iter!=pToSort->forwardNeighborsEnd(i);iter++) {
            if(iter->second<0){
                numberOfNegative++;
            }

        }
    }
    size_t numberOfPositive=pToSort->getNumberOfEdges()-numberOfNegative;

    size_t positiveSetSize=size_t(std::round(parameters.getPositiveThresholdLifted()*numberOfPositive));
    size_t negativeSetSize=size_t(std::round(parameters.getNegativeThresholdLifted()*numberOfNegative));

    //TODO if zero -> return zero
    std::multiset<double> positiveCosts;
    std::multiset<double> negativeCosts;
    for (int i = 0; i < pToSort->getNumberOfVertices(); ++i) {
        auto iter=pToSort->forwardNeighborsBegin(i);
        for (;iter!=pToSort->forwardNeighborsEnd(i);iter++) {
            double cost=iter->second;

            if(cost<0){
                if(negativeCosts.size()<negativeSetSize){
                    negativeCosts.insert(cost);
                }
                else{
                    auto iter2=negativeCosts.begin();
                    if(*iter2<cost){
                        //  std::cout<<"neg costs size before erase "<<negativeCosts.size()<<std::endl;
                        negativeCosts.erase(iter2);
                          //std::cout<<"neg costs size after erase "<<negativeCosts.size()<<std::endl;
                        negativeCosts.insert(cost);
                         // std::cout<<"neg costs size after insert "<<negativeCosts.size()<<std::endl;
                    }
                    //std::cout<<"neg costs size "<<negativeCosts.size()<<std::endl;

                }
            }
            else{
                if(positiveCosts.size()<positiveSetSize){
                    positiveCosts.insert(cost);
                }
                else{
                    auto iter2=positiveCosts.end();
                    iter--;
                    if(*iter2>cost){
                        positiveCosts.erase(iter2);
                        positiveCosts.insert(cost);
                    }

                }
            }

        }
    }

    assert(negativeCosts.size()==negativeSetSize);
    assert(positiveCosts.size()==positiveSetSize);
    positiveLiftedThreshold=*positiveCosts.rbegin();
    negativeLiftedThreshold=*negativeCosts.begin();

    pToSort=pBaseGraph;
    size_t baseSetSize=size_t(std::round(parameters.getBaseUpperThreshold()*pBaseGraph->getNumberOfEdges()));

    std::set<double> baseCosts;
    for (int i = 0; i < pToSort->getNumberOfVertices(); ++i) {
        auto iter=pToSort->forwardNeighborsBegin(i);
        for (;iter!=pToSort->forwardNeighborsEnd(i);iter++) {
            double cost=iter->second;

            if(baseCosts.size()<baseSetSize){
                baseCosts.insert(cost);
            }
            else{
                auto iter=baseCosts.begin();
                if(*iter<cost){
                    baseCosts.erase(iter);
                    baseCosts.insert(cost);
                }

            }
        }
    }

    baseThreshold=*baseCosts.begin();

  if(diagnostics())std::cout<<"negative threshold lifted "<<negativeLiftedThreshold<<std::endl;

    if(diagnostics())std::cout<<"positive threshold lifted "<<positiveLiftedThreshold<<std::endl;


    if(diagnostics())std::cout<<"base thresholds "<<baseThreshold<<", index "<<std::endl;

}





LdpInstance::LdpInstance(LdpParameters<>& configParameters,const py::array_t<size_t>& baseEdges,const py::array_t<size_t>& liftedEdges,const  py::array_t<double>& baseCosts,const  py::array_t<double>& liftedCosts,const  py::array_t<double>& verticesCosts,VertexGroups<>& pvg):parameters(configParameters){

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    pCompleteGraph=nullptr;

    const auto baseEdgeVector=baseEdges.unchecked<2>();
    const std::size_t dimBase1=baseEdgeVector.shape(0);
    const std::size_t dimBase2=baseEdgeVector.shape(1);
    const auto baseCostVector=baseCosts.unchecked<1>();
    const size_t dimOfBaseCosts=baseCostVector.shape(0);

    const auto liftedEdgeVector=liftedEdges.unchecked<2>();
    const std::size_t dimLifted1=liftedEdgeVector.shape(0);
    const std::size_t dimLifted2=liftedEdgeVector.shape(1);
    const auto liftedCostVector=liftedCosts.unchecked<1>();
    const size_t dimOfLiftedCosts=liftedCostVector.shape(0);



    if(dimBase2!=2||dimLifted2!=2){
        std::string message="Wrong dimension of edge array, second dimension 2 expected";
        throw std::invalid_argument(message);
    }
    if(dimBase1!=dimOfBaseCosts||dimLifted1!=dimOfLiftedCosts){
        std::string message="Dimension of edge array and edge costs do not match.";
        throw std::invalid_argument(message);
    }


    vertexGroups=pvg;
    myGraph=LdpDirectedGraph(baseEdgeVector,baseCostVector,parameters.getInputCost(),parameters.getOutputCost());
    myGraphLifted=LdpDirectedGraph(liftedEdgeVector,liftedCostVector);
    s_=myGraphLifted.getNumberOfEdges();
    t_=s_+1;


    if(configParameters.isSparsify()){
        if(parameters.isUseAdaptiveThreshold()){
           //initAdaptiveThresholds(&csBase.completeScore,&csLifted.completeScore);
            initAdaptiveThresholds(&myGraph,&myGraphLifted);
        }
        else{
            baseThreshold=parameters.getBaseUpperThreshold();
            positiveLiftedThreshold=parameters.getPositiveThresholdLifted();
            negativeLiftedThreshold=parameters.getNegativeThresholdLifted();
        }
        sparsifyBaseGraphNew(myGraph);
        reachable=initReachableLdp(myGraph,parameters,&vertexGroups);
        //disjointPaths::keepFractionOfLifted(*this,configParameters);
        sparsifyLiftedGraphNew(myGraphLifted);
    }
    else{
        reachable=initReachableLdp(myGraph,parameters,&vertexGroups);
    }



    initLiftedStructure();

    numberOfVertices=myGraph.getNumberOfVertices();
    numberOfLiftedEdges=myGraphLifted.getNumberOfEdges();
    numberOfEdges=myGraph.getNumberOfEdges();

   // vertexScore=std::vector<double>(graph_.numberOfVertices()-2);
    const auto costVector=verticesCosts.unchecked<1>();
    const size_t dimOfCosts=costVector.shape(0);
    vertexScore=std::vector<double>(myGraphLifted.getNumberOfVertices());
    assert(dimOfCosts==vertexScore.size());
    for (size_t i=0;i<dimOfCosts;i++){
        vertexScore[i]=costVector(i);
    }



    assert(vertexScore.size()==numberOfVertices-2);

    minV=0;
    maxV=vertexScore.size();  //TODO: Check usage!

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    parameters.getControlOutput() << "Time of instance constructor = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << " seconds" << std::endl;
    parameters.writeControlOutput();

}

void LdpInstance::init(){


   if(diagnostics())  std::cout<<"number of vertices "<<myGraph.getNumberOfVertices()<<std::endl;
   strongBaseEdges=std::vector<std::unordered_set<size_t>>(myGraph.getNumberOfVertices());
    if(parameters.isSparsify()){


        sparsifyBaseGraphNew(myGraph,parameters.isAllBaseZero());

          reachable=initReachableLdp(myGraph,parameters,&vertexGroups);

        sparsifyLiftedGraphNew(myGraphLifted);


        initLiftedStructure();


    }
    else{
        std::cout<<"Initialization of base and lifted graph from one graph without sparsification not supported"<<std::endl;
        assert(false);
        reachable=initReachableLdp(myGraph,parameters,&vertexGroups);
        initLiftedStructure();
    }

    numberOfEdges=myGraph.getNumberOfEdges();
    numberOfLiftedEdges=myGraphLifted.getNumberOfEdges();

}




double LdpInstance::evaluateClustering(const std::vector<size_t>& labels) const{

    double value=0;
    size_t maxLabel=0;
    if(pCompleteGraph!=nullptr){
        assert(pCompleteGraph->getNumberOfVertices()==labels.size());
        for (int i = 0; i < labels.size(); ++i) {
            if(labels[i]!=0){
                maxLabel=std::max(maxLabel,labels[i]);
                auto iter=pCompleteGraph->forwardNeighborsBegin(i);
                for (;iter!=pCompleteGraph->forwardNeighborsEnd(i);iter++) {
                    size_t vertex=iter->first;
                    double cost=iter->second;
                    assert(vertex<labels.size());
                    if(labels[vertex]==labels[i]) value+=cost;
                }
            }

        }
    }
    if(parameters.isMustCutMissing()){
        size_t mustCuts=0;
        std::vector<std::vector<size_t>> verticesToLabels(maxLabel);
        for (size_t i = 0; i < labels.size(); ++i) {
            if(labels[i]>0){
                size_t labelOrder=labels[i]-1;
                assert(labelOrder<verticesToLabels.size());
                size_t size=verticesToLabels[labelOrder].size();
                for (size_t j = 0; j < size; ++j) {
                    size_t vertex=verticesToLabels[labelOrder][j];
                    if(isReachable(vertex,i)&&!canJoin(vertex,i)){
                        mustCuts++;
                    }
                }

            }

        }
        if(diagnostics()&&mustCuts>0){
            std::cout<<"number of must cuts: "<<mustCuts<<std::endl;
        }
        value+=100*mustCuts;
    }
    return value;

}





void LdpInstance::initLiftedStructure(){
    size_t n=numberOfVertices-2;

    //std::cout<<"number of vertices "<<numberOfVertices<<std::endl;
    sncNeighborStructure=std::vector<size_t>(n+2);
    sncBUNeighborStructure=std::vector<size_t>(n+2);
    sncTDStructure=std::vector<double>(n+2);
    sncBUStructure=std::vector<double>(n+2);
    sncClosedVertices=std::vector<char>(n+2);
    sncLiftedMessages=std::vector<double>(n+2);
    sncVerticesInScope=std::vector<char>(n+2);

    liftedStructure=std::vector<ShiftedVector<char>>(n);
    for (size_t i = 0; i < n; ++i) {
        size_t maxVertex=0;
        size_t minVertex=n;
        auto iter=myGraphLifted.forwardNeighborsBegin(i);
        for(;iter!=myGraphLifted.forwardNeighborsEnd(i);iter++){
            size_t vertex2=iter->first;

           if(vertex2<minVertex){
               minVertex=vertex2;
           }
           if(vertex2>maxVertex){
               maxVertex=vertex2;
           }
        }
        liftedStructure[i]=ShiftedVector<char>(minVertex,maxVertex,false);

        iter=myGraphLifted.forwardNeighborsBegin(i);
        for(;iter!=myGraphLifted.forwardNeighborsEnd(i);iter++){
            size_t vertex2=iter->first;
           liftedStructure[i][vertex2]=1;
        }
    }
}




void LdpInstance::initCanJoinStructure(const LdpDirectedGraph& completeGraph){
    size_t n=completeGraph.getNumberOfVertices();


    canJoinStructure=std::vector<ShiftedVector<char>>(n);
    for (size_t i = 0; i < n; ++i) {
        size_t maxVertex=0;
        size_t minVertex=n;
        auto iter=completeGraph.forwardNeighborsBegin(i);
        for(;iter!=completeGraph.forwardNeighborsEnd(i);iter++){
            size_t vertex2=iter->first;

           if(vertex2<minVertex){
               minVertex=vertex2;
           }
           if(vertex2>maxVertex){
               maxVertex=vertex2;
           }
        }
        canJoinStructure[i]=ShiftedVector<char>(minVertex,maxVertex,false);

        iter=completeGraph.forwardNeighborsBegin(i);
        for(;iter!=completeGraph.forwardNeighborsEnd(i);iter++){
            size_t vertex2=iter->first;
            canJoinStructure[i][vertex2]=1;
        }
    }
}







//void LdpInstance::sparsifyLiftedGraphNew(const LdpDirectedGraph& inputLiftedGraph){


//    parameters.getControlOutput()<<"Sparsify lifted graph"<<std::endl;
//    parameters.writeControlOutput();
//    //TODO run automaticLifted to find candidates first

//    assert(parameters.getMaxTimeLifted()>=1);

//    std::vector<std::array<size_t,2>> edges;
//    std::vector<double> costs;

//    for (size_t i = 0; i < inputLiftedGraph.getNumberOfVertices(); ++i) {
//        size_t l0=vertexGroups.getGroupIndex(i);
//        auto iter=inputLiftedGraph.forwardNeighborsBegin(i);
//        for (;iter!=inputLiftedGraph.forwardNeighborsEnd(i);iter++) {
//            size_t w=iter->first;
//            size_t l1=vertexGroups.getGroupIndex(w);
//            assert(l1>l0);
//            if(isReachable(i,w)){
//                double cost=iter->second;
//                if(cost<negativeLiftedThreshold||cost>positiveLiftedThreshold){ //Previously was strictly higher or lower
//                    if(l1-l0==1){
//                        if(!parameters.isAllBaseZero()){
//                            if(parameters.getDenseTimeLifted()>=1||parameters.getLongerIntervalLifted()==1){
//                                std::pair<size_t,double>* baseIt=myGraph.forwardNeighborsBegin(i);
//                                while(baseIt->first!=w){
//                                    baseIt++;
//                                    assert(baseIt!=myGraph.forwardNeighborsEnd(i));
//                                }
//                                //std::cout<<"base edge with cost "<<i<<" "<<w<<std::endl;
//                                baseIt->second+=cost;
//                                baseIt=myGraph.backwardNeighborsBegin(w);
//                                while(baseIt->first!=i){
//                                    baseIt++;
//                                    assert(baseIt!=myGraph.forwardNeighborsEnd(w));
//                                }
//                                //std::cout<<"base edge with cost "<<i<<" "<<w<<std::endl;
//                                baseIt->second+=cost;
//                                //std::cout<<"base edge with cost "<<i<<" "<<w<<": "<<baseIt->second<<std::endl;
//                            }
//                        }
//                        //TODO do not use lifted, add cost to base if max time lifted is geq 1
//                    }
//                    else{
//                        bool useEdge=l1-l0<=parameters.getDenseTimeLifted();
//                        if(!useEdge){
//                            size_t timeGapDiff=l1-l0-parameters.getDenseTimeLifted();
//                            useEdge=((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0);
//                        }
//                        if(useEdge){
//                            //std::cout<<"lifted edge "<<i<<" "<<w<<": "<<cost<<std::endl;
//                            edges.push_back({i,w});
//                            costs.push_back(cost);
//                        }


//                    }
//                }
//            }
//        }
//    }

//    EdgeVector ev(edges);
//    InfoVector iv(costs);
//    myGraphLifted=LdpDirectedGraph(ev,iv);

//    parameters.getControlOutput()<<"left "<<myGraphLifted.getNumberOfEdges()<<" lifted edges"<<std::endl;


//    parameters.writeControlOutput();
//    initLiftedStructure();

//}



void LdpInstance::sparsifyLiftedGraphNew(const LdpDirectedGraph& inputLiftedGraph){
//TODO remove input parameter

    parameters.getControlOutput()<<"Sparsify lifted graph"<<std::endl;
    parameters.writeControlOutput();
    //TODO run automaticLifted to find candidates first

    assert(parameters.getMaxTimeLifted()>=1);

    std::vector<std::array<size_t,2>> edges;
    std::vector<double> costs;

    size_t highCostEdges=0;
    size_t maxTimeGapInGraph=1;

    for (size_t i = 0; i < inputLiftedGraph.getNumberOfVertices(); ++i) {
        size_t l0=vertexGroups.getGroupIndex(i);
        auto iterLifted=inputLiftedGraph.forwardNeighborsBegin(i);
        auto iterBase=myGraph.forwardNeighborsBegin(i);
        auto baseEnd=myGraph.forwardNeighborsEnd(i);
        //        for (;iterLifted!=inputLiftedGraph.forwardNeighborsEnd(i);iterLifted++) {
        while(iterLifted!=inputLiftedGraph.forwardNeighborsEnd(i)) {
            while(iterBase!=baseEnd&&iterBase->first<iterLifted->first){
                iterBase++;
            }
            bool isSame=iterBase!=baseEnd&&iterLifted->first==iterBase->first;
            size_t w=iterLifted->first;

            size_t l1=vertexGroups.getGroupIndex(w);
            // if(w>=1356) std::cout<<"candidate edge "<<i<<","<<w<<" cost "<<iterLifted->second<<", gap "<<(l1-l0)<<std::endl;
            maxTimeGapInGraph=std::max(maxTimeGapInGraph,size_t(l1-l0));
            assert(l1>l0);
            if(isReachable(i,w)){
            //  std::cout<<"reachable"<<std::endl;
                double cost=iterLifted->second;
                //if(cost<negativeLiftedThreshold||cost>positiveLiftedThreshold){ //Previously was strictly higher or lower
                if(l1-l0==1){
                    if(cost<negativeLiftedThreshold||cost>positiveLiftedThreshold){
                        if(!parameters.isAllBaseZero()){
                            if(isSame&&(parameters.getDenseTimeLifted()>=1||parameters.getLongerIntervalLifted()==1)){
                                iterBase->second+=cost;
                                std::pair<size_t,double>*baseIt=myGraph.backwardNeighborsBegin(w);
                                while(baseIt->first!=i){
                                    baseIt++;
                                    assert(baseIt!=myGraph.forwardNeighborsEnd(w));
                                }
                                //std::cout<<"base edge with cost "<<i<<" "<<w<<std::endl;
                                baseIt->second+=cost;
                                //std::cout<<"base edge with cost "<<i<<" "<<w<<": "<<baseIt->second<<std::endl;
                            }
                        }
                    }

                }
                else{
                   // bool useEdge=parameters.isAllBaseZero()&&isSame;
                    //if(!useEdge){
                    bool useEdge=false;
                        bool goodCost=(cost<negativeLiftedThreshold||cost>positiveLiftedThreshold);
                        if(goodCost){
                            // if(w>=1356) std::cout<<"has good cost"<<std::endl;
                            useEdge=(l1-l0<=parameters.getDenseTimeLifted());

                            if(!useEdge){
                                size_t timeGapDiff=l1-l0-parameters.getDenseTimeLifted();
                                useEdge=((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0);
                            }
                        }
                        if(useEdge){
                            //if(w>=1356) std::cout<<"lifted edge "<<i<<" "<<w<<": "<<cost<<std::endl;
                            edges.push_back({i,w});
                            costs.push_back(cost);
                            if(cost>=100.00){
                                highCostEdges++;
                            }

                        }
                        else if(isSame&&parameters.isAllBaseZero()){
                            if(parameters.isBaseCoverdWithLifted()){ //add lifted edge
                                edges.push_back({i,w});
                                costs.push_back(cost);
                               // if(w>=1356) std::cout<<"edge "<<i<<","<<w<<std::endl;
                            }
                            else{  //add base cost
                                iterBase->second+=cost;
                                std::pair<size_t,double>*baseIt=myGraph.backwardNeighborsBegin(w);
                                while(baseIt->first!=i){
                                    baseIt++;
                                    assert(baseIt!=myGraph.forwardNeighborsEnd(w));
                                }
                                //std::cout<<"base edge with cost "<<i<<" "<<w<<std::endl;
                                baseIt->second+=cost;
                            }
                        }


                }

            }
            else{
              //  if(w>=1356) std::cout<<"not reachable"<<std::endl;
            }
            iterLifted++;
        }
    }

    parameters.getControlOutput()<<"number of high cost edges "<<highCostEdges<<std::endl;

    parameters.writeControlOutput();

    size_t reachableMustCut=0;
    maxTimeGapInGraph=std::min(maxTimeGapInGraph,parameters.getMaxTimeLifted());

    if(parameters.isMustCutMissing()){
        for (size_t i = 0; i < inputLiftedGraph.getNumberOfVertices(); ++i) {
            for(auto& d:reachable[i]){
                if(d!=i&&d<inputLiftedGraph.getNumberOfVertices()&&!canJoin(i,d)){
                    size_t l0=vertexGroups.getGroupIndex(i);
                    size_t l1=vertexGroups.getGroupIndex(d);
                    if((l1-l0)>1&&(l1-l0)<=maxTimeGapInGraph){
                        if((l1-l0)<=parameters.getDenseTimeLifted()||((l1-l0-parameters.getDenseTimeLifted())%parameters.getLongerIntervalLifted()==0)){
                            edges.push_back({i,d});
                            costs.push_back(parameters.getMustCutPenalty());
                            reachableMustCut++;
                            //std::cout<<i<<", "<<d<<", time diff "<<(l1-l0)<<std::endl;
                        }
                    }
                }

            }

        }
    }
   parameters.getControlOutput()<<"used reachable must cuts: "<<reachableMustCut<<std::endl;

    EdgeVector ev(edges);
    InfoVector iv(costs);
    myGraphLifted=LdpDirectedGraph(ev,iv,inputLiftedGraph.getNumberOfVertices());

    parameters.getControlOutput()<<"Lifted vertices "<<myGraphLifted.getNumberOfVertices()<<". Left "<<myGraphLifted.getNumberOfEdges()<<" lifted edges"<<std::endl;


    parameters.writeControlOutput();
    initLiftedStructure();

}

void LdpInstance::sparsifyBaseGraphNew(const LdpDirectedGraph& inputGraph, bool zeroCost){
    parameters.getControlOutput()<<"Sparsify base graph"<<std::endl;
    parameters.writeControlOutput();

    std::vector<double> newBaseCosts;

    size_t k=parameters.getKnnK();

    std::vector<std::array<size_t,2>> edgesToUse;
    std::vector<double> costsToUse;

    size_t n=inputGraph.getNumberOfVertices()-2;
    //std::cout<<"base sparse n "<<n<<std::endl;

    //std::vector<std::vector<std::multimap<double,size_t>>> edgesToKeep(n);
    std::vector<std::vector<std::set<size_t>>> edgesToKeep(n);
    std::vector<std::vector<std::multimap<double,size_t>>> edgesToKeepBackward(n);

    for (size_t i = 0; i < n; ++i) {
        edgesToKeep[i]=std::vector<std::set<size_t>>(parameters.getMaxTimeBase());
        edgesToKeepBackward[i]=std::vector<std::multimap<double,size_t>>(parameters.getMaxTimeBase());
    }

    for (size_t v0 = 0; v0 < inputGraph.getNumberOfVertices(); ++v0) {
        //std::cout<<"vertex 0 "<<v0<<std::endl;
        std::vector<std::multimap<double,size_t>> edgesToKeepForward(parameters.getMaxTimeBase());
        size_t l0=vertexGroups.getGroupIndex(v0);
        //std::cout<<"layer "<<l0<<std::endl;
        auto iter=inputGraph.forwardNeighborsBegin(v0);
        for (; iter!=inputGraph.forwardNeighborsEnd(v0); iter++) {

            size_t v1=iter->first;

            double cost=iter->second;

            if(v0==s_||v1==t_){

                edgesToUse.push_back({v0,v1});


                costsToUse.push_back(cost);
            }
            else{
                size_t l1=vertexGroups.getGroupIndex(v1);
                assert(l1>l0);
                size_t gap=l1-l0; //map to vector index
                //assert(gap<=parameters.getMaxTimeBase());

                if(gap<=parameters.getMaxTimeBase()){
                    if(gap<=parameters.getKnnTimeGap()||cost<=parameters.getBaseUpperThreshold()){
                        assert(edgesToKeep.size()>gap-1);
                        if(edgesToKeepForward[gap-1].size()<parameters.getKnnK()){
                            edgesToKeepForward[gap-1].insert(std::pair<double,size_t>(cost,v1));
                        }
                        else{
                            double lastValue=edgesToKeepForward[gap-1].rbegin()->first;
                            if(cost<lastValue){
                                auto it=edgesToKeepForward[gap-1].end();
                                it--;
                                edgesToKeepForward[gap-1].erase(it);
                                edgesToKeepForward[gap-1].insert(std::pair<double,size_t>(cost,v1));
                            }
                        }

                        if(edgesToKeepBackward[v1][gap-1].size()<parameters.getKnnK()){
                            edgesToKeepBackward[v1][gap-1].insert(std::pair<double,size_t>(cost,v0));
                        }
                        else{
                            double lastValue=edgesToKeepBackward[v1][gap-1].rbegin()->first;
                            if(cost<lastValue){
                                auto it=edgesToKeepBackward[v1][gap-1].end();
                                it--;
                                edgesToKeepBackward[v1][gap-1].erase(it);
                                edgesToKeepBackward[v1][gap-1].insert(std::pair<double,size_t>(cost,v0));
                            }
                        }

                    }
                }
            }
        }
        for (size_t i = 0; i < parameters.getMaxTimeBase(); ++i) {

            for (auto it=edgesToKeepForward[i].begin();it!=edgesToKeepForward[i].end();it++) {
                edgesToKeep[v0][i].insert(it->second);
            }

        }
    }

     for (size_t v1 = 0; v1 < n; ++v1) {
        for (size_t i = 0; i < edgesToKeepBackward[v1].size(); ++i) {
            auto it=edgesToKeepBackward[v1][i].rbegin();
            for (;it!=edgesToKeepBackward[v1][i].rend();it++) {
                double cost=it->first;
                size_t v0=it->second;

                if(edgesToKeep[v0][i].count(v1)>0){

                    edgesToUse.push_back({v0,v1});
                    if(!zeroCost||i==0){
                        costsToUse.push_back(cost);
                    }
                    else{
                        costsToUse.push_back(0.0);
                    }
                }

            }

        }

     }

    EdgeVector ev(edgesToUse);
    InfoVector iv(costsToUse);
    myGraph=LdpDirectedGraph(ev,iv);

    parameters.getControlOutput()<<"Base vertices "<<myGraph.getNumberOfVertices()<<". Left "<<myGraph.getNumberOfEdges()<<" base edges."<<std::endl;
    parameters.writeControlOutput();

 }




bool LdpInstance::checkStrongBase(const size_t &v, const size_t &w) const{
    auto iter=myGraph.forwardNeighborsBegin(v);
    bool alternativePath=false;
    while(!alternativePath&&iter!=myGraph.forwardNeighborsEnd(v)&&iter->first<w){
        if(isReachable(iter->first,w)){
            alternativePath=true;
        }
        iter++;
    }
    return !alternativePath;
}





bool LdpInstance::isStrongBase(size_t v,size_t w) const{
    bool isStrong= strongBaseEdges.at(v).count(w)>0;
    return isStrong;
   // return false;
}




std::vector<std::unordered_set<size_t>> LdpInstance::initReachableLdp(const LdpDirectedGraph & graph,LdpParameters<>& parameters,const VertexGroups<size_t>* vg){
    parameters.getControlOutput()<<"Run Floyd Warshall"<<std::endl;
    const size_t n=graph.getNumberOfVertices();


    std::vector<std::unordered_set<size_t>> desc(n);
    std::vector<std::vector<std::bitset<10000>>> descBit(n);
    size_t columns=n/10000;
    size_t mod=n%10000;
    if(mod>0) columns++;


    for (size_t i = 0; i < n; ++i) {
        descBit[i]=std::vector<std::bitset<10000>>(columns);
    }

    for (size_t v = 0; v < n; ++v) {

        descBit[v][v/10000][v%10000]=1; //make this reflexive
        auto iter=graph.forwardNeighborsBegin(v);
        for (;iter!=graph.forwardNeighborsEnd(v);iter++) {
            size_t w=iter->first;
            descBit[v][w/10000][w%10000]=1;
        }

    }



    if(vg==nullptr){
        for (size_t k1 = 0; k1 <columns; ++k1) {
            for (size_t k2 = 0; k2 < 10000; ++k2) {
                if(k1*10000+k2<n){
                    for (size_t i = 0; i < n; ++i) {
                        if(descBit[i][k1][k2]){
                            for (size_t j = 0; j < columns; ++j) {
                                descBit[i][j]|=descBit[k1*10000+k2][j];
                            }

                        }
                    }
                }
                else{
                    break;
                }
            }

        }
    }
    else{
        for (size_t k1 = 0; k1 <columns; ++k1) {
            for (size_t k2 = 0; k2 < 10000; ++k2) {
                if(k1*10000+k2<n){
                    size_t maxTime=vg->getGroupIndex(k1*10000+k2);
                    size_t minTime=0;
                    if(maxTime>parameters.getMaxTimeGapComplete()){
                        minTime=maxTime-parameters.getMaxTimeGapComplete();
                    }

                    for (size_t t = minTime; t < maxTime; ++t) {//TODO use time gap
                        const std::vector<size_t>& vertices=vg->getGroupVertices(t);
                        for (size_t i:vertices) {
                            assert(i<descBit.size());
                            if(descBit[i][k1][k2]){
                                for (size_t j = 0; j < columns; ++j) {
                                    descBit[i][j]|=descBit[k1*10000+k2][j];
                                }

                            }
                        }
                    }
                }
                else{
                    break;
                }
            }

        }

    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t k1 = 0; k1 <columns; ++k1) {
            for (size_t k2 = 0; k2 < 10000; ++k2) {
                if(k1*10000+k2<n&&descBit[i][k1][k2]){
                    desc[i].insert(k1*10000+k2);
                }
            }
        }

    }

    size_t reachableMustCut=0;
    if(parameters.isMustCutMissing()){
        for (size_t i = 0; i < canJoinStructure.size(); ++i) {
            for(auto& d:desc[i]){
                if(vg!=nullptr&&d<canJoinStructure.size()){
                    size_t l0=vg->getGroupIndex(i);
                    size_t l1=vg->getGroupIndex(d);
                    if(l1!=l0&&(l1-l0)<=parameters.getMaxTimeGapComplete()&&!canJoin(i,d)){
                        reachableMustCut++;
                    }
                }
            }

        }
    }
   parameters.getControlOutput()<<"number of reachable must cuts: "<<reachableMustCut<<std::endl;

   parameters.writeControlOutput();
    return desc;

}



}}//End of namespaces
