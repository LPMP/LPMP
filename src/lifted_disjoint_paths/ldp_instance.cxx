#include"lifted_disjoint_paths/ldp_instance.hxx"


namespace LPMP{
namespace lifted_disjoint_paths {






LdpInstance::LdpInstance(LdpParameters<> &configParameters, CompleteStructure<>& cs):
    parameters(configParameters)
{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


    //size_t minT=1;
    //size_t maxT=cs.maxTime+1;
    //readGraphWithTime(minT,maxT,&cs);

    std::cout<<"instance constructor "<<std::endl;
    myGraphLifted=cs.myCompleteGraph;
    std::cout<<"complete graph assigned "<<myGraphLifted.getNumberOfVertices()<<std::endl;
    s_=myGraphLifted.getNumberOfVertices();
    t_=s_+1;
    myGraph=LdpDirectedGraph(cs.myCompleteGraph,parameters.getInputCost(),parameters.getOutputCost());


    std::cout<<"graphs init"<<std::endl;
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
//    if(debug()) std::cout<<"Adding automatic lifted edges"<<std::endl;
//    for (size_t i = 0; i < graph_.numberOfEdges(); ++i) {
//        size_t v0=graph_.vertexOfEdge(i,0);
//        size_t v1=graph_.vertexOfEdge(i,1);
//        if(v0!=s_&&v1!=t_){
//            //	if(secOrderDesc[v0][v1]){
//            graphLifted_.insertEdge(v0,v1);
//            liftedEdgeScore.push_back(edgeScore[i]);

//        }
//    }
    numberOfVertices=myGraph.getNumberOfVertices();

    std::cout<<"before method init"<<std::endl;

    init();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const std::chrono::steady_clock::time_point& csBegin=cs.getContructorBegin();
    parameters.getControlOutput()<< "Time of instance constructor = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << " seconds" << std::endl;
    parameters.getControlOutput()<< "Time of instance constructor and complete structure = " << std::chrono::duration_cast<std::chrono::seconds> (end - csBegin).count() << " seconds" << std::endl;
    parameters.writeControlOutput();




}

LdpInstance::LdpInstance(LdpParameters<>& configParameters,LdpBatchProcess& BP):
    parameters(configParameters){

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //std::cout<<"instance constructor"<<std::endl;
    //graph_=BP.getOutputGraph();
    //edgeScore=BP.getEdgeScore();
    vertexScore=BP.getVerticesScore();
    //std::cout<<"edges and graph set, creating vg"<<std::endl;
    BP.createLocalVG(vertexGroups);
    //std::cout<<"created vg"<<std::endl;

    myGraphLifted=BP.getMyCompleteGraph();
    s_=myGraphLifted.getNumberOfVertices();
    t_=s_+1;
    myGraph=LdpDirectedGraph(BP.getMyCompleteGraph(),parameters.getInputCost(),parameters.getOutputCost());


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
   // if(diagnostics()) std::cout << "Time of instance constructor = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << " seconds" << std::endl;
     parameters.getControlOutput() << "Time of instance constructor = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << " seconds" << std::endl;
    const std::chrono::steady_clock::time_point& bpBegin=BP.getContructorBegin();
    parameters.getControlOutput() << "Time of instance constructor and batch process = " << std::chrono::duration_cast<std::chrono::seconds> (end - bpBegin).count() << " seconds" << std::endl;
    parameters.writeControlOutput();


}






//LdpInstance::LdpInstance(LdpParameters<>& configParameters,LdpBatchProcess& BP):
//    parameters(configParameters){

//    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//    //std::cout<<"instance constructor"<<std::endl;
//    //graph_=BP.getOutputGraph();
//    edgeScore=BP.getEdgeScore();
//    vertexScore=BP.getVerticesScore();
//    //std::cout<<"edges and graph set, creating vg"<<std::endl;
//    BP.createLocalVG(vertexGroups);
//    //std::cout<<"created vg"<<std::endl;


//    numberOfVertices=BP.getOutputGraph().numberOfVertices()+2;
//    graph_=andres::graph::Digraph<>(numberOfVertices);
//    for (size_t i = 0; i < BP.getOutputGraph().numberOfEdges(); ++i) {
//        graph_.insertEdge(BP.getOutputGraph().vertexOfEdge(i,0),BP.getOutputGraph().vertexOfEdge(i,1));
//    }
//    s_=numberOfVertices-2;
//    t_=s_+1;


//    for (size_t v = 0; v < numberOfVertices-2; ++v) {
//        graph_.insertEdge(s_,v);
//        edgeScore.push_back(parameters.getInputCost());
//        graph_.insertEdge(v,t_);
//        edgeScore.push_back(parameters.getOutputCost());
//    }

//    graphLifted_ = BP.getOutputGraph();
//    liftedEdgeScore=BP.getEdgeScore();

//    if(parameters.isUseAdaptiveThreshold()){
//        if(diagnostics()) std::cout<<"using adaptive"<<std::endl;
//       initAdaptiveThresholds(&edgeScore,nullptr);
//    }
//    else{
//        baseThreshold=parameters.getBaseUpperThreshold();
//        positiveLiftedThreshold=parameters.getPositiveThresholdLifted();
//        negativeLiftedThreshold=parameters.getNegativeThresholdLifted();
//    }


//    init();

//    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//   // if(diagnostics()) std::cout << "Time of instance constructor = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << " seconds" << std::endl;
//     parameters.getControlOutput() << "Time of instance constructor = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << " seconds" << std::endl;
//    const std::chrono::steady_clock::time_point& bpBegin=BP.getContructorBegin();
//    parameters.getControlOutput() << "Time of instance constructor and batch process = " << std::chrono::duration_cast<std::chrono::seconds> (end - bpBegin).count() << " seconds" << std::endl;
//    parameters.writeControlOutput();


//}


//LdpInstance::LdpInstance( LdpParameters<>& configParameters):
//    parameters(configParameters)
//{
//    char delim=',';
//    size_t maxVertex;

//        vertexGroups=disjointPaths::VertexGroups<size_t>(parameters);
//        maxVertex=vertexGroups.getMaxVertex();

//    std::ifstream graphFile(parameters.getGraphFileName());
//    readGraph(graphFile,maxVertex,delim);
//    graphFile.close();

//    init();

//    //baseEdgeLabels=std::vector<bool>(numberOfEdges);
//}

//void LdpInstance::initAdaptiveThresholds(const std::vector<double>* baseCosts,const std::vector<double>* liftedCosts){
//    std::vector<double> costsToSort;
//    if(liftedCosts!=nullptr){
//        costsToSort=*liftedCosts;
//    }
//    else{
//        costsToSort=*baseCosts;
//    }
//    assert(parameters.getPositiveThresholdLifted()>=0&&parameters.getPositiveThresholdLifted()<1);
//    assert(parameters.getNegativeThresholdLifted()>=0&&parameters.getNegativeThresholdLifted()<1);
//    assert(parameters.getBaseUpperThreshold()>=0&&parameters.getBaseUpperThreshold()<1);
//    std::sort(costsToSort.begin(),costsToSort.end());
//    size_t numberOfNegative=0;
//    while(costsToSort[numberOfNegative]<0){
//        numberOfNegative++;
//    }
//    if(diagnostics()) std::cout<<"number of negative "<<numberOfNegative<<std::endl;
//    size_t numberOfPositive=costsToSort.size()-numberOfNegative;

//    if(diagnostics()) std::cout<<"number of positive "<<numberOfPositive<<std::endl;


//    size_t indexForNegThreshold=size_t(std::round(numberOfNegative*(1.0-parameters.getNegativeThresholdLifted())));
//    if(indexForNegThreshold>0) indexForNegThreshold--;
//    assert(indexForNegThreshold<costsToSort.size());
//    negativeLiftedThreshold=costsToSort[indexForNegThreshold];
//    if(diagnostics())std::cout<<"negative threshold lifted "<<negativeLiftedThreshold<<", index "<<indexForNegThreshold<<std::endl;

//    size_t indexForPositiveThreshold=size_t(std::round(numberOfPositive*(parameters.getPositiveThresholdLifted())))+numberOfNegative;
//    if(indexForPositiveThreshold>0) indexForPositiveThreshold--;
//    assert(indexForPositiveThreshold<costsToSort.size());
//    positiveLiftedThreshold=costsToSort[indexForPositiveThreshold];

//    if(diagnostics())std::cout<<"positive threshold lifted "<<positiveLiftedThreshold<<", index "<<indexForPositiveThreshold<<std::endl;

//    if(liftedCosts!=nullptr){
//        costsToSort=*baseCosts;
//        std::sort(costsToSort.begin(),costsToSort.end());
//    }
//    size_t indexForBaseThreshold=size_t(std::round(costsToSort.size()*(1.0-parameters.getBaseUpperThreshold())));
//    if(indexForBaseThreshold>0) indexForBaseThreshold--;
//    assert(indexForBaseThreshold<costsToSort.size());
//    baseThreshold=costsToSort[indexForBaseThreshold];

//    if(diagnostics())std::cout<<"base thresholds "<<baseThreshold<<", index "<<indexForBaseThreshold<<std::endl;

//}


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
    std::set<double> positiveCosts;
    std::set<double> negativeCosts;
    for (int i = 0; i < pToSort->getNumberOfVertices(); ++i) {
        auto iter=pToSort->forwardNeighborsBegin(i);
        for (;iter!=pToSort->forwardNeighborsEnd(i);iter++) {
            double cost=iter->second;

            if(cost<0){
                if(negativeCosts.size()<negativeSetSize){
                    negativeCosts.insert(cost);
                }
                else{
                    auto iter=negativeCosts.begin();
                    if(*iter<cost){
                        negativeCosts.erase(iter);
                        negativeCosts.insert(cost);
                    }

                }
            }
            else{
                if(positiveCosts.size()<positiveSetSize){
                    positiveCosts.insert(cost);
                }
                else{
                    auto iter=positiveCosts.end();
                    iter--;
                    if(*iter>cost){
                        positiveCosts.erase(iter);
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



//LdpInstance::LdpInstance(LdpParameters<>& configParameters,const std::vector<std::array<size_t,2>>& completeEdges,const  std::vector<double>& completeCosts,disjointPaths::VertexGroups<>& pvg)
//    :parameters(configParameters)
//{

//    vertexGroups=pvg;
//    assert(completeCosts.size()==completeEdges.size());
//    std::vector<double> costsToSort=completeCosts;

//    if(configParameters.isUseAdaptiveThreshold()){


//    }



//}


//LdpInstance::LdpInstance(LdpParameters<>& configParameters,const disjointPaths::TwoGraphsInputStructure& twoGraphsIS):parameters(configParameters){
LdpInstance::LdpInstance(LdpParameters<>& configParameters,const py::array_t<size_t>& baseEdges,const py::array_t<size_t>& liftedEdges,const  py::array_t<double>& baseCosts,const  py::array_t<double>& liftedCosts,const  py::array_t<double>& verticesCosts,VertexGroups<>& pvg):parameters(configParameters){
//LdpInstance::LdpInstance(LdpParameters<>& configParameters,const std::vector<std::array<size_t,2>>& baseEdges,const std::vector<std::array<size_t,2>>& liftedEdges,const  std::vector<double>& baseCosts,const  std::vector<double>& liftedCosts,disjointPaths::VertexGroups<>& pvg):parameters(configParameters){
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


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

//    numberOfVertices=myGraph.getNumberOfVertices();
//    assert(numberOfVertices==t_+1);

//    csBase.addEdgesFromVectorsAll(baseEdges,baseCosts);
//    size_t maxVertex=vertexGroups.getMaxVertex();
//    graph_=andres::graph::Digraph<>(maxVertex+3);
//   // graph_.reserveEdges(csBase.completeGraph.numberOfEdges()+2*csBase.completeGraph.numberOfVertices());
//    for (size_t i = 0; i < csBase.completeGraph.numberOfEdges(); ++i) {
//         size_t v=csBase.completeGraph.vertexOfEdge(i,0);
//         size_t w=csBase.completeGraph.vertexOfEdge(i,1);
//         graph_.insertEdge(v,w);
//    }
//    edgeScore=csBase.completeScore;


//    s_=maxVertex+1;
//    t_=s_+1;
//    for (size_t i = 0; i <= maxVertex; ++i) {
//        graph_.insertEdge(s_,i);
//        edgeScore.push_back(configParameters.getInputCost());
//        graph_.insertEdge(i,t_);
//        edgeScore.push_back(configParameters.getOutputCost());
//    }

//    numberOfVertices=graph_.numberOfVertices();

//    CompleteStructure<> csLifted(vertexGroups);
//    csLifted.setVerticesCosts(verticesCosts);
//    csLifted.addEdgesFromVectorsAll(liftedEdges,liftedCosts);

//    graphLifted_=csLifted.completeGraph;
//    liftedEdgeScore=csLifted.completeScore;


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

//    myGraph=LdpDirectedGraph(graph_,edgeScore);
//    myGraphLifted=LdpDirectedGraph(graphLifted_,liftedEdgeScore);

    minV=0;
    maxV=vertexScore.size();  //TODO: Check usage!

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    parameters.getControlOutput() << "Time of instance constructor = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << " seconds" << std::endl;
    parameters.writeControlOutput();

}

void LdpInstance::init(){


    if(debug()) std::cout<<"done"<<std::endl;

   if(diagnostics())  std::cout<<"number of vertices "<<myGraph.getNumberOfVertices()<<std::endl;
   strongBaseEdges=std::vector<std::unordered_set<size_t>>(myGraph.getNumberOfVertices());
    if(parameters.isSparsify()){


        std::cout<<"before sparsify"<<std::endl;
        sparsifyBaseGraphNew(myGraph,parameters.isAllBaseZero());
      //  reachable=initReachableSet(graph_,parameters,&vertexGroups);
          reachable=initReachableLdp(myGraph,parameters,&vertexGroups);

//        if(parameters.isAllBaseZero()){
//            //std::cout<<"base to zero"<<std::endl;
//            std::fill(edgeScore.begin(),edgeScore.end(),0);
//        }
//        else{
//            std::cout<<"base cost preserved"<<std::endl;
//        }
        sparsifyLiftedGraphNew(myGraphLifted);

        std::cout<<"after sparsify"<<std::endl;
        initLiftedStructure();


    }
    else{
        std::cout<<"Initialization of base and lifted graph from one graph without sparsification not supported"<<std::endl;
        assert(false);
        reachable=initReachableLdp(myGraph,parameters,&vertexGroups);
        //reachable=initReachableSet(graph_,parameters,&vertexGroups);
        initLiftedStructure();
    }

//    myGraph=LdpDirectedGraph(graph_,edgeScore);
//    myGraphLifted=LdpDirectedGraph(graphLifted_,liftedEdgeScore);

//    numberOfVertices=graph_.numberOfVertices();
//    numberOfEdges=graph_.numberOfEdges();
//    numberOfLiftedEdges=graphLifted_.numberOfEdges();
    numberOfEdges=myGraph.getNumberOfEdges();
    numberOfLiftedEdges=myGraphLifted.getNumberOfEdges();

}



//void LdpInstance::readGraphWithTime(size_t minTime,size_t maxTime,CompleteStructure<>* cs){

//	andres::graph::Digraph<>& completeGraph=cs->completeGraph;
//	std::vector<double>& completeScore=cs->completeScore;
//    VertexGroups<> vg=cs->getVertexGroups();

//	std::unordered_map<size_t,std::vector<size_t>> groups;

//	size_t mt=minTime;
//	while(vg.getGroupVertices(mt).size()==0){
//		mt++;
//	}
//	size_t minVertex=vg.getGroupVertices(mt)[0];

//	mt=maxTime-1;
//	while(vg.getGroupVertices(mt).size()==0){
//		mt--;
//	}
//	size_t maxVertex=*(vg.getGroupVertices(mt).rbegin());

//	size_t numberOfVertices=maxVertex-minVertex+3;
//	s_ = numberOfVertices - 2;
//	t_ = numberOfVertices - 1;

//	std::vector<size_t> vToGroup(numberOfVertices);
//	vToGroup[s_]=0;
//	vToGroup[t_]=maxTime-minTime+1;

//    vertexScore = std::vector<double>(numberOfVertices, 0);
//    const std::vector<double>& verticesScoreComplete=cs->verticesScore;
//	for (int gi = minTime; gi < maxTime; ++gi) {
//		//groups[gi-minTime+1]=std::vector<size_t>();
//		for(size_t v:vg.getGroupVertices(gi)){
//			size_t vertex=v-minVertex;
//            assert(v<verticesScoreComplete.size());
//            assert(vertex<vertexScore.size());
//            vertexScore[vertex]=verticesScoreComplete[v];
//			groups[gi-minTime+1].push_back(vertex);
//			vToGroup[vertex]=gi-minTime+1;
//		}
//	}

//    vertexGroups=VertexGroups<>(groups,vToGroup);

//	graphLifted_ = andres::graph::Digraph<>(numberOfVertices);
//	graph_ = andres::graph::Digraph<>(numberOfVertices);

//	bool useZeroInOut=false;
//	for (int v = 0; v < numberOfVertices-2; ++v) {
//		graph_.insertEdge(s_,v);
//		if(useZeroInOut) edgeScore.push_back(0);
//		else edgeScore.push_back(parameters.getInputCost());
//		graph_.insertEdge(v,t_);
//		if(useZeroInOut) edgeScore.push_back(0);
//		else edgeScore.push_back(parameters.getOutputCost());
//	}


//	for (int v = minVertex; v < maxVertex; ++v) {
//		for (int i = 0; i < completeGraph.numberOfEdgesFromVertex(v); ++i) {
//			size_t w=completeGraph.vertexFromVertex(v,i);
//			if(w>maxVertex) continue;
//			size_t e=completeGraph.edgeFromVertex(v,i);
//			graph_.insertEdge(v-minVertex, w-minVertex);
//			edgeScore.push_back(completeScore[e]);
//		}
//	}
//	minV=minVertex;
//	maxV=maxVertex;

//}





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
//        for (size_t j = 0; j < graphLifted_.numberOfEdgesFromVertex(i); ++j) {
//           size_t vertex2=graphLifted_.vertexFromVertex(i,j);
           if(vertex2<minVertex){
               minVertex=vertex2;
           }
           if(vertex2>maxVertex){
               maxVertex=vertex2;
           }
        }
        liftedStructure[i]=ShiftedVector<char>(minVertex,maxVertex,false);
        //for (size_t j = 0; j < graphLifted_.numberOfEdgesFromVertex(i); ++j) {
         //  size_t vertex2=graphLifted_.vertexFromVertex(i,j);
        iter=myGraphLifted.forwardNeighborsBegin(i);
        for(;iter!=myGraphLifted.forwardNeighborsEnd(i);iter++){
            size_t vertex2=iter->first;
           liftedStructure[i][vertex2]=1;
        }
    }
}




//void LdpInstance::sparsifyLiftedGraph(){


//    parameters.getControlOutput()<<"Sparsify lifted graph"<<std::endl;
//    parameters.writeControlOutput();
//    //TODO run automaticLifted to find candidates first

//    std::vector<double> newLiftedCosts;

//    andres::graph::Digraph<> tempGraphLifted=(graph_.numberOfVertices());


//    std::unordered_map<size_t,std::set<size_t>> liftedEdges;
//    for (size_t v = 0; v < graphLifted_.numberOfVertices()-2; ++v) {
//        std::unordered_set<size_t> alternativePath;
//        for (size_t i = 0; i < graph_.numberOfEdgesFromVertex(v); ++i) {
//            size_t w=graph_.vertexFromVertex(v,i);
//            for(size_t u:reachable[w]){
//                if(u!=w){
//                    alternativePath.insert(u);
//                }
//            }
//        }
//        for (size_t i = 0; i < graph_.numberOfEdgesFromVertex(v); ++i) {
//            size_t w=graph_.vertexFromVertex(v,i);
//            if(alternativePath.count(w)==0){
//                strongBaseEdges.at(v).insert(w);
//            }
//        }

//        for (size_t i = 0; i < graphLifted_.numberOfEdgesFromVertex(v); ++i) {
//            size_t w=graphLifted_.vertexFromVertex(v,i);
//            if(w!=t_){
//                if(alternativePath.count(w)>0) liftedEdges[v].insert(w);

//            }
//        }
//    }



//    parameters.getControlOutput()<<"done"<<std::endl;
//   parameters.writeControlOutput();


//    for (int i = 0; i < graphLifted_.numberOfEdges(); ++i) {
//        size_t v0=graphLifted_.vertexOfEdge(i,0);
//        size_t v1=graphLifted_.vertexOfEdge(i,1);
//        int l0=vertexGroups.getGroupIndex(v0);
//        int l1=vertexGroups.getGroupIndex(v1);
//        double cost=getLiftedEdgeScore(i);
//        bool goodCost=(cost<negativeLiftedThreshold)||(cost>positiveLiftedThreshold);
//        if(isReachable(v0,v1)){

//            int timeGapDiff=l1-l0-parameters.getDenseTimeLifted();
//            bool timeConstraint=l1-l0<=parameters.getDenseTimeLifted()||((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0);
//            if(timeConstraint&&goodCost){
//                if(liftedEdges[v0].count(v1)>0){
//                    tempGraphLifted.insertEdge(v0,v1);
//                    newLiftedCosts.push_back(cost);
//                }
//                else{
//                    auto edgeTest=graph_.findEdge(v0,v1);
//                    if(edgeTest.first){
//                        edgeScore[edgeTest.second]+=cost;  //Compensate that the lifted edge has been removed
//                    }

//                }

//            }
//        }

//    }



//    liftedEdgeScore=newLiftedCosts;

//    graphLifted_=tempGraphLifted;
//      parameters.getControlOutput()<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
//    parameters.writeControlOutput();
//    if(graphLifted_.numberOfEdges()!=newLiftedCosts.size()){
//         parameters.getControlOutput()<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
//        parameters.getControlOutput()<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
//    }
//    else{
//         parameters.getControlOutput()<<"lifted edge number and lifted graph size match "<<std::endl;
//        parameters.getControlOutput()<<"lifted edge number and lifted graph size match "<<std::endl;

//    }
//   parameters.writeControlOutput();
//    initLiftedStructure();

//}



void LdpInstance::sparsifyLiftedGraphNew(const LdpDirectedGraph& inputLiftedGraph){


    parameters.getControlOutput()<<"Sparsify lifted graph"<<std::endl;
    parameters.writeControlOutput();
    //TODO run automaticLifted to find candidates first

    assert(parameters.getMaxTimeLifted()>=1);

    std::vector<std::array<size_t,2>> edges;
    std::vector<double> costs;

    for (size_t i = 0; i < inputLiftedGraph.getNumberOfVertices(); ++i) {
        size_t l0=vertexGroups.getGroupIndex(i);
        auto iter=inputLiftedGraph.forwardNeighborsBegin(i);
        for (;iter!=inputLiftedGraph.forwardNeighborsEnd(i);iter++) {
            size_t w=iter->first;
            size_t l1=vertexGroups.getGroupIndex(w);
            assert(l1>l0);
            if(isReachable(i,w)){
                double cost=iter->second;
                if(cost<=negativeLiftedThreshold||cost>=positiveLiftedThreshold){ //Previously was strictly higher or lower
                    if(l1-l0==1){
                        if(parameters.getDenseTimeLifted()>=1||parameters.getLongerIntervalLifted()==1){
                            std::pair<size_t,double>* baseIt=myGraph.forwardNeighborsBegin(i);
                            while(baseIt->first!=w){
                                baseIt++;
                                assert(baseIt!=myGraph.forwardNeighborsEnd(i));
                            }
                            baseIt->second+=cost;
                        }
                        //TODO do not use lifted, add cost to base if max time lifted is geq 1
                    }
                    else{
                        bool useEdge=l1-l0<=parameters.getDenseTimeLifted();
                        if(!useEdge){
                            size_t timeGapDiff=l1-l0-parameters.getDenseTimeLifted();
                            useEdge=((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0);
                        }
                        if(useEdge){
                            edges.push_back({i,w});
                            costs.push_back(cost);
                        }


                    }
                }
            }
        }
    }

    EdgeVector ev(edges);
    InfoVector iv(costs);
    myGraphLifted=LdpDirectedGraph(ev,iv);

    parameters.getControlOutput()<<"left "<<myGraphLifted.getNumberOfEdges()<<" lifted edges"<<std::endl;


    parameters.writeControlOutput();
    initLiftedStructure();

}








////It could get complete structure and directly output two dim array
//void LdpInstance::sparsifyBaseGraphNew(andres::graph::Digraph<> &inputGraph){
//     parameters.getControlOutput()<<"Sparsify base graph"<<std::endl;
//     parameters.writeControlOutput();

//     std::vector<double> newBaseCosts;
//     //std::vector<size_t> inOutEdges;
//     size_t k=parameters.getKnnK();
//     //std::vector<size_t> goodLongEdges;

//   //  andres::graph::Digraph<> tempGraph(graph_.numberOfVertices());


//     size_t nrVertices=inputGraph.numberOfVertices();
//     std::vector<std::array<size_t,2>> edgesForMyGraph;

//     std::vector<bool> finalEdges(inputGraph.numberOfEdges(),false);
//     for (int v0 = 0; v0 < nrVertices; ++v0) {
//         std::unordered_map<int,std::list<size_t>> edgesToKeep;
//         size_t l0=vertexGroups.getGroupIndex(v0);
//         for (size_t ne = 0; ne < inputGraph.numberOfEdgesFromVertex(v0); ++ne) {
//             size_t e=inputGraph.edgeFromVertex(v0,ne);
//             size_t v1=inputGraph.vertexFromVertex(v0,ne);
// //			std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<problemGraph.getEdgeScore(e)<<std::endl;
//             if(v0==s_||v1==t_){
//                 //tempGraph.insertEdge(v0,v1);
//                 //newBaseCosts.push_back(edgeScore[e]);
//                 finalEdges[e]=true;
//             }
//             else{
//                 size_t l1=vertexGroups.getGroupIndex(v1);
//                 size_t gap=l1-l0;
//                 if(gap<=parameters.getMaxTimeBase()){
//                 //if(gap<=parameters.getKnnTimeGap()){
//                     //gap=std::min(parameters.getKnnTimeGap()+1,gap);
//                     double cost=edgeScore[e];
//                     if(edgesToKeep.count(gap)>0){
//                         std::list<size_t>& smallList=edgesToKeep[gap];
//                         auto it=smallList.begin();
//                         double bsf=edgeScore[*it];
//                         //std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
//                         while(bsf>cost&&it!=smallList.end()){
//                             it++;
//                             size_t index=*it;
//                             if(it!=smallList.end()){
//                                 bsf=edgeScore[index];
//                                 //	std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
//                             }
//                         }
//                         if(it!=smallList.begin()){
//                             smallList.insert(it,e);
//                             if(smallList.size()>k) smallList.pop_front();
//                         }
//                         else if(smallList.size()<k){
//                             smallList.push_front(e);
//                         }
//                     }
//                     else{
//                         edgesToKeep[gap].push_front(e);
//                     }
//                 }
// //				else if(gap<=parameters.getMaxTimeBase()){
// //					if(getEdgeScore(e)<=parameters.getBaseUpperThreshold()){
// //						//tempGraph.insertEdge(v0,v1);
// //						//newBaseCosts.push_back(getEdgeScore(e));
// //						finalEdges[e]=true;
// //					}
// //
// //				}
//             }
//         }
//         //std::cout.precision(4);
//         double bsf=0;
//         for (int gap = 0; gap <= parameters.getKnnTimeGap(); ++gap) {
//             if(edgesToKeep.count(gap)>0){
//                 auto& smallList=edgesToKeep[gap];
//                 for(size_t e:smallList){
//                     finalEdges[e]=true;
//                     if(edgeScore[e]<bsf){
//                         bsf=edgeScore[e];
//                     }
//                 }
//             }
//         }

//         for (int gap =  parameters.getKnnTimeGap()+1;gap<=parameters.getMaxTimeBase(); ++gap) {
//             if(edgesToKeep.count(gap)>0){
//                 auto& smallList=edgesToKeep[gap];
//                 for(size_t e:smallList){
//                     double score=edgeScore[e];
//                     if(score<=baseThreshold){
//                         finalEdges[e]=true;
//                     }
//                 }

//             }
//         }

//     }

//     for (int e = 0; e < inputGraph.numberOfEdges(); ++e) {
//         if(finalEdges[e]){
//             size_t v0=inputGraph.vertexOfEdge(e,0);
//             size_t v1=inputGraph.vertexOfEdge(e,1);
//             //tempGraph.insertEdge(v0,v1);
//             edgesForMyGraph.push_back({v0,v1});
//             newBaseCosts.push_back(edgeScore[e]);
//         }
//     }

//     assert(edgesForMyGraph.size()==newBaseCosts.size());

//     edgeScore=newBaseCosts;


//     if(edgesForMyGraph.size()!=newBaseCosts.size()){
//         parameters.getControlOutput()<<"edge number mismatch, edge vector: "<<edgesForMyGraph.size()<<", cost vector "<<newBaseCosts.size()<<std::endl;
//         parameters.writeControlOutput();

//     }
//     else{

//         parameters.getControlOutput()<<"edge number and graph size match "<<std::endl;
//         parameters.writeControlOutput();
//     }

//     parameters.getControlOutput()<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
//     parameters.writeControlOutput();

//     EdgeVector ev(edgesForMyGraph);
//     InfoVector iv(newBaseCosts);
//     //myGraph=LdpDirectedGraph(edgesForMyGraph,newBaseCosts);
//     myGraph=LdpDirectedGraph(ev,iv);

// }



//It could get complete structure and directly output two dim array
//void LdpInstance::sparsifyBaseGraph(){
//     parameters.getControlOutput()<<"Sparsify base graph"<<std::endl;
//     parameters.writeControlOutput();

//     std::vector<double> newBaseCosts;
//     //std::vector<size_t> inOutEdges;
//     size_t k=parameters.getKnnK();
//     //std::vector<size_t> goodLongEdges;

//     andres::graph::Digraph<> tempGraph(graph_.numberOfVertices());



//     std::vector<bool> finalEdges(graph_.numberOfEdges(),false);
//     for (int v0 = 0; v0 < graph_.numberOfVertices(); ++v0) {
//         std::unordered_map<int,std::list<size_t>> edgesToKeep;
//         size_t l0=vertexGroups.getGroupIndex(v0);
//         for (size_t ne = 0; ne < graph_.numberOfEdgesFromVertex(v0); ++ne) {
//             size_t e=graph_.edgeFromVertex(v0,ne);
//             size_t v1=graph_.vertexFromVertex(v0,ne);
// //			std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<problemGraph.getEdgeScore(e)<<std::endl;
//             if(v0==s_||v1==t_){
//                 //tempGraph.insertEdge(v0,v1);
//                 //newBaseCosts.push_back(edgeScore[e]);
//                 finalEdges[e]=true;
//             }
//             else{
//                 size_t l1=vertexGroups.getGroupIndex(v1);
//                 size_t gap=l1-l0;
//                 if(gap<=parameters.getMaxTimeBase()){
//                 //if(gap<=parameters.getKnnTimeGap()){
//                     //gap=std::min(parameters.getKnnTimeGap()+1,gap);
//                     double cost=edgeScore[e];
//                     if(edgesToKeep.count(gap)>0){
//                         std::list<size_t>& smallList=edgesToKeep[gap];
//                         auto it=smallList.begin();
//                         double bsf=edgeScore[*it];
//                         //std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
//                         while(bsf>cost&&it!=smallList.end()){
//                             it++;
//                             size_t index=*it;
//                             if(it!=smallList.end()){
//                                 bsf=edgeScore[index];
//                                 //	std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
//                             }
//                         }
//                         if(it!=smallList.begin()){
//                             smallList.insert(it,e);
//                             if(smallList.size()>k) smallList.pop_front();
//                         }
//                         else if(smallList.size()<k){
//                             smallList.push_front(e);
//                         }
//                     }
//                     else{
//                         edgesToKeep[gap].push_front(e);
//                     }
//                 }
// //				else if(gap<=parameters.getMaxTimeBase()){
// //					if(getEdgeScore(e)<=parameters.getBaseUpperThreshold()){
// //						//tempGraph.insertEdge(v0,v1);
// //						//newBaseCosts.push_back(getEdgeScore(e));
// //						finalEdges[e]=true;
// //					}
// //
// //				}
//             }
//         }
//         //std::cout.precision(4);
//         double bsf=0;
//         for (int gap = 0; gap <= parameters.getKnnTimeGap(); ++gap) {
//             if(edgesToKeep.count(gap)>0){
//                 auto& smallList=edgesToKeep[gap];
//                 for(size_t e:smallList){
//                     finalEdges[e]=true;
//                     if(edgeScore[e]<bsf){
//                         bsf=edgeScore[e];
//                     }
//                 }
//             }
//         }

//         for (int gap =  parameters.getKnnTimeGap()+1;gap<=parameters.getMaxTimeBase(); ++gap) {
//             if(edgesToKeep.count(gap)>0){
//                 auto& smallList=edgesToKeep[gap];
//                 for(size_t e:smallList){
//                     double score=edgeScore[e];
//                     if(score<=baseThreshold){
//                         finalEdges[e]=true;
//                     }
//                 }

//             }
//         }

//     }

//     for (int e = 0; e < graph_.numberOfEdges(); ++e) {
//         if(finalEdges[e]){
//             size_t v0=graph_.vertexOfEdge(e,0);
//             size_t v1=graph_.vertexOfEdge(e,1);
//             tempGraph.insertEdge(v0,v1);
//             newBaseCosts.push_back(edgeScore[e]);
//         }
//     }

//     if(newBaseCosts.size()!=tempGraph.numberOfEdges()){
//         throw std::runtime_error("Error in base graph sparsification.");
//     }



//     graph_=tempGraph;
//     edgeScore=newBaseCosts;


//     if(graph_.numberOfEdges()!=newBaseCosts.size()){
//         parameters.getControlOutput()<<"edge number mismatch, graph: "<<graph_.numberOfEdges()<<", cost vector "<<newBaseCosts.size()<<std::endl;
//         parameters.writeControlOutput();

//     }
//     else{

//         parameters.getControlOutput()<<"edge number and graph size match "<<std::endl;
//         parameters.writeControlOutput();
//     }

//     parameters.getControlOutput()<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
//     parameters.writeControlOutput();

// }



void LdpInstance::sparsifyBaseGraphNew(const LdpDirectedGraph& inputGraph, bool zeroCost){
    parameters.getControlOutput()<<"Sparsify base graph"<<std::endl;
    parameters.writeControlOutput();

    std::vector<double> newBaseCosts;
    //std::vector<size_t> inOutEdges;
    size_t k=parameters.getKnnK();
    //std::vector<size_t> goodLongEdges;

  //  andres::graph::Digraph<> tempGraph(graph_.numberOfVertices());


    std::vector<std::array<size_t,2>> edgesToUse;
    std::vector<double> costsToUse;


    for (size_t v0 = 0; v0 < inputGraph.getNumberOfVertices(); ++v0) {
        std::vector<std::multimap<double,size_t>> edgesToKeep(parameters.getMaxTimeBase());
        size_t l0=vertexGroups.getGroupIndex(v0);
        auto iter=inputGraph.forwardNeighborsBegin(v0);
        for (; iter!=inputGraph.forwardNeighborsEnd(v0); iter++) {
            size_t v1=iter->first;
            double cost=iter->second;
            //			std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<problemGraph.getEdgeScore(e)<<std::endl;
            if(v0==s_||v1==t_){
                //tempGraph.insertEdge(v0,v1);
                //newBaseCosts.push_back(edgeScore[e]);
                edgesToUse.push_back({v0,v1});
                costsToUse.push_back(cost);
            }
            else{
                size_t l1=vertexGroups.getGroupIndex(v1);
                assert(l1>l0);
                size_t gap=l1-l0; //map to vector index
                assert(gap<parameters.getMaxTimeBase());
                if(gap<=parameters.getKnnTimeGap()||cost<=parameters.getBaseUpperThreshold()){
                    if(edgesToKeep[gap-1].size()<parameters.getKnnK()){
                        edgesToKeep[gap-1].insert(std::pair<double,size_t>(cost,v1));
                    }
                    else{
                        double lastValue=edgesToKeep[gap-1].rbegin()->second;
                        if(cost<lastValue){
                            auto it=edgesToKeep[gap-1].end();
                            it--;
                            edgesToKeep[gap-1].erase(it);
                            edgesToKeep[gap-1].insert(std::pair<double,size_t>(cost,v1));
                        }
                    }

                }
            }
        }
        for (size_t i = 0; i < edgesToKeep.size(); ++i) {
            auto it=edgesToKeep[i].begin();
            for (;it!=edgesToKeep[i].end();it++) {
                double cost=it->first;
                size_t v1=it->second;
                edgesToUse.push_back({v0,v1});
                if(zeroCost){
                    costsToUse.push_back(0.0);
                }
                else{
                    costsToUse.push_back(cost);
                }
            }

        }

     }

    EdgeVector ev(edgesToUse);
    InfoVector iv(costsToUse);
    myGraph=LdpDirectedGraph(ev,iv);

 }









bool LdpInstance::isStrongBase(size_t v,size_t w) const{
    bool isStrong= strongBaseEdges.at(v).count(w)>0;
    return isStrong;
   // return false;
}




std::vector<std::unordered_set<size_t>> LdpInstance::initReachableLdp(const LdpDirectedGraph & graph,LdpParameters<>& parameters,const VertexGroups<size_t>* vg){


    //        levinkov::Timer tfw;
    //                        tfw.start();


    parameters.getControlOutput()<<"Run Floyd Warshall"<<std::endl;
    const size_t n=graph.getNumberOfVertices();

    // todo: use some matrix structure
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

    return desc;

}



}}//End of namespaces
