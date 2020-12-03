#include"lifted_disjoint_paths/ldp_instance.hxx"


namespace LPMP{
namespace lifted_disjoint_paths {






LdpInstance::LdpInstance(LdpParameters<> &configParameters, CompleteStructure<>& cs):
    parameters(configParameters)
{

    size_t minT=1;
    size_t maxT=cs.maxTime+1;
    readGraphWithTime(minT,maxT,&cs);


    if(parameters.isUseAdaptiveThreshold()){
        if(diagnostics()) std::cout<<"using adaptive"<<std::endl;
       initAdaptiveThresholds(&cs.completeScore,nullptr);
    }
    else{
        baseThreshold=parameters.getBaseUpperThreshold();
        positiveLiftedThreshold=parameters.getPositiveThresholdLifted();
        negativeLiftedThreshold=parameters.getNegativeThresholdLifted();
    }
    if(debug()) std::cout<<"Adding automatic lifted edges"<<std::endl;
    for (size_t i = 0; i < graph_.numberOfEdges(); ++i) {
        size_t v0=graph_.vertexOfEdge(i,0);
        size_t v1=graph_.vertexOfEdge(i,1);
        if(v0!=s_&&v1!=t_){
            //	if(secOrderDesc[v0][v1]){
            graphLifted_.insertEdge(v0,v1);
            liftedEdgeScore.push_back(edgeScore[i]);

        }
    }
    numberOfVertices=graph_.numberOfVertices();

    init();

}


LdpInstance::LdpInstance(LdpParameters<>& configParameters,LdpBatchProcess& BP):
    parameters(configParameters){
    //std::cout<<"instance constructor"<<std::endl;
    //graph_=BP.getOutputGraph();
    edgeScore=BP.getEdgeScore();
    vertexScore=BP.getVerticesScore();
    //std::cout<<"edges and graph set, creating vg"<<std::endl;
    BP.createLocalVG(vertexGroups);
    //std::cout<<"created vg"<<std::endl;


    numberOfVertices=BP.getOutputGraph().numberOfVertices()+2;
    graph_=andres::graph::Digraph<>(numberOfVertices);
    for (size_t i = 0; i < BP.getOutputGraph().numberOfEdges(); ++i) {
        graph_.insertEdge(BP.getOutputGraph().vertexOfEdge(i,0),BP.getOutputGraph().vertexOfEdge(i,1));
    }
    s_=numberOfVertices-2;
    t_=s_+1;


    for (size_t v = 0; v < numberOfVertices-2; ++v) {
        graph_.insertEdge(s_,v);
        edgeScore.push_back(parameters.getInputCost());
        graph_.insertEdge(v,t_);
        edgeScore.push_back(parameters.getOutputCost());
    }

    graphLifted_ = BP.getOutputGraph();
    liftedEdgeScore=BP.getEdgeScore();

    if(parameters.isUseAdaptiveThreshold()){
        if(diagnostics()) std::cout<<"using adaptive"<<std::endl;
       initAdaptiveThresholds(&edgeScore,nullptr);
    }
    else{
        baseThreshold=parameters.getBaseUpperThreshold();
        positiveLiftedThreshold=parameters.getPositiveThresholdLifted();
        negativeLiftedThreshold=parameters.getNegativeThresholdLifted();
    }


    init();

}


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

void LdpInstance::initAdaptiveThresholds(const std::vector<double>* baseCosts,const std::vector<double>* liftedCosts){
    std::vector<double> costsToSort;
    if(liftedCosts!=nullptr){
        costsToSort=*liftedCosts;
    }
    else{
        costsToSort=*baseCosts;
    }
    assert(parameters.getPositiveThresholdLifted()>=0&&parameters.getPositiveThresholdLifted()<1);
    assert(parameters.getNegativeThresholdLifted()>=0&&parameters.getNegativeThresholdLifted()<1);
    assert(parameters.getBaseUpperThreshold()>=0&&parameters.getBaseUpperThreshold()<1);
    std::sort(costsToSort.begin(),costsToSort.end());
    size_t numberOfNegative=0;
    while(costsToSort[numberOfNegative]<0){
        numberOfNegative++;
    }
    if(diagnostics()) std::cout<<"number of negative "<<numberOfNegative<<std::endl;
    size_t numberOfPositive=costsToSort.size()-numberOfNegative;

    if(diagnostics()) std::cout<<"number of positive "<<numberOfPositive<<std::endl;


    size_t indexForNegThreshold=size_t(std::round(numberOfNegative*(1.0-parameters.getNegativeThresholdLifted())));
    if(indexForNegThreshold>0) indexForNegThreshold--;
    assert(indexForNegThreshold<costsToSort.size());
    negativeLiftedThreshold=costsToSort[indexForNegThreshold];
    if(diagnostics())std::cout<<"negative threshold lifted "<<negativeLiftedThreshold<<", index "<<indexForNegThreshold<<std::endl;

    size_t indexForPositiveThreshold=size_t(std::round(numberOfPositive*(parameters.getPositiveThresholdLifted())))+numberOfNegative;
    if(indexForPositiveThreshold>0) indexForPositiveThreshold--;
    assert(indexForPositiveThreshold<costsToSort.size());
    positiveLiftedThreshold=costsToSort[indexForPositiveThreshold];

    if(diagnostics())std::cout<<"positive threshold lifted "<<positiveLiftedThreshold<<", index "<<indexForPositiveThreshold<<std::endl;

    if(liftedCosts!=nullptr){
        costsToSort=*baseCosts;
        std::sort(costsToSort.begin(),costsToSort.end());
    }
    size_t indexForBaseThreshold=size_t(std::round(costsToSort.size()*(1.0-parameters.getBaseUpperThreshold())));
    if(indexForBaseThreshold>0) indexForBaseThreshold--;
    assert(indexForBaseThreshold<costsToSort.size());
    baseThreshold=costsToSort[indexForBaseThreshold];

    if(diagnostics())std::cout<<"base thresholds "<<baseThreshold<<", index "<<indexForBaseThreshold<<std::endl;

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
    CompleteStructure<> csBase(pvg);

    vertexGroups=pvg;


    csBase.addEdgesFromVectorsAll(baseEdges,baseCosts);
    size_t maxVertex=vertexGroups.getMaxVertex();
    graph_=andres::graph::Digraph<>(maxVertex+3);
   // graph_.reserveEdges(csBase.completeGraph.numberOfEdges()+2*csBase.completeGraph.numberOfVertices());
    for (size_t i = 0; i < csBase.completeGraph.numberOfEdges(); ++i) {
         size_t v=csBase.completeGraph.vertexOfEdge(i,0);
         size_t w=csBase.completeGraph.vertexOfEdge(i,1);
         graph_.insertEdge(v,w);
    }
    edgeScore=csBase.completeScore;


    s_=maxVertex+1;
    t_=s_+1;
    for (size_t i = 0; i <= maxVertex; ++i) {
        graph_.insertEdge(s_,i);
        edgeScore.push_back(configParameters.getInputCost());
        graph_.insertEdge(i,t_);
        edgeScore.push_back(configParameters.getOutputCost());
    }

    numberOfVertices=graph_.numberOfVertices();

    CompleteStructure<> csLifted(vertexGroups);
    csLifted.setVerticesCosts(verticesCosts);
    csLifted.addEdgesFromVectorsAll(liftedEdges,liftedCosts);

    graphLifted_=csLifted.completeGraph;
    liftedEdgeScore=csLifted.completeScore;


    if(configParameters.isSparsify()){
        if(parameters.isUseAdaptiveThreshold()){
           initAdaptiveThresholds(&csBase.completeScore,&csLifted.completeScore);
        }
        else{
            baseThreshold=parameters.getBaseUpperThreshold();
            positiveLiftedThreshold=parameters.getPositiveThresholdLifted();
            negativeLiftedThreshold=parameters.getNegativeThresholdLifted();
        }
        sparsifyBaseGraph();
        reachable=initReachableSet(graph_,parameters,&vertexGroups);
        //disjointPaths::keepFractionOfLifted(*this,configParameters);
        sparsifyLiftedGraph();
    }
    else{
        reachable=initReachableSet(graph_,parameters,&vertexGroups);
    }



    initLiftedStructure();

    numberOfVertices=graph_.numberOfVertices();
    numberOfLiftedEdges=graphLifted_.numberOfEdges();
    numberOfEdges=graph_.numberOfEdges();

   // vertexScore=std::vector<double>(graph_.numberOfVertices()-2);
    vertexScore=csLifted.verticesScore;
    assert(vertexScore.size()==numberOfVertices-2);

    myGraph=LdpDirectedGraph(graph_,edgeScore);
    myGraphLifted=LdpDirectedGraph(graphLifted_,liftedEdgeScore);

    minV=0;
    maxV=maxVertex;


}

void LdpInstance::init(){



    if(debug()) std::cout<<"done"<<std::endl;

   if(diagnostics())  std::cout<<"number of vertices "<<graph_.numberOfVertices()<<std::endl;
   strongBaseEdges=std::vector<std::unordered_set<size_t>>(graph_.numberOfVertices());
    if(parameters.isSparsify()){


        sparsifyBaseGraph();
        reachable=initReachableSet(graph_,parameters,&vertexGroups);

        if(parameters.isAllBaseZero()){
            //std::cout<<"base to zero"<<std::endl;
            std::fill(edgeScore.begin(),edgeScore.end(),0);
        }
//        else{
//            std::cout<<"base cost preserved"<<std::endl;
//        }
        sparsifyLiftedGraph();
        initLiftedStructure();


    }
    else{
        std::cout<<"Initialization of base and lifted graph from one graph without sparsification not supported"<<std::endl;
        assert(false);

        reachable=initReachableSet(graph_,parameters,&vertexGroups);
        initLiftedStructure();
    }

    myGraph=LdpDirectedGraph(graph_,edgeScore);
    myGraphLifted=LdpDirectedGraph(graphLifted_,liftedEdgeScore);

    numberOfVertices=graph_.numberOfVertices();
    numberOfEdges=graph_.numberOfEdges();
    numberOfLiftedEdges=graphLifted_.numberOfEdges();

}



void LdpInstance::readGraphWithTime(size_t minTime,size_t maxTime,CompleteStructure<>* cs){

	andres::graph::Digraph<>& completeGraph=cs->completeGraph;
	std::vector<double>& completeScore=cs->completeScore;
    VertexGroups<> vg=cs->getVertexGroups();

	std::unordered_map<size_t,std::vector<size_t>> groups;

	size_t mt=minTime;
	while(vg.getGroupVertices(mt).size()==0){
		mt++;
	}
	size_t minVertex=vg.getGroupVertices(mt)[0];

	mt=maxTime-1;
	while(vg.getGroupVertices(mt).size()==0){
		mt--;
	}
	size_t maxVertex=*(vg.getGroupVertices(mt).rbegin());

	size_t numberOfVertices=maxVertex-minVertex+3;
	s_ = numberOfVertices - 2;
	t_ = numberOfVertices - 1;

	std::vector<size_t> vToGroup(numberOfVertices);
	vToGroup[s_]=0;
	vToGroup[t_]=maxTime-minTime+1;

    vertexScore = std::vector<double>(numberOfVertices, 0);
    const std::vector<double>& verticesScoreComplete=cs->verticesScore;
	for (int gi = minTime; gi < maxTime; ++gi) {
		//groups[gi-minTime+1]=std::vector<size_t>();
		for(size_t v:vg.getGroupVertices(gi)){
			size_t vertex=v-minVertex;
            assert(v<verticesScoreComplete.size());
            assert(vertex<vertexScore.size());
            vertexScore[vertex]=verticesScoreComplete[v];
			groups[gi-minTime+1].push_back(vertex);
			vToGroup[vertex]=gi-minTime+1;
		}
	}

    vertexGroups=VertexGroups<>(groups,vToGroup);

	graphLifted_ = andres::graph::Digraph<>(numberOfVertices);
	graph_ = andres::graph::Digraph<>(numberOfVertices);

	bool useZeroInOut=false;
	for (int v = 0; v < numberOfVertices-2; ++v) {
		graph_.insertEdge(s_,v);
		if(useZeroInOut) edgeScore.push_back(0);
		else edgeScore.push_back(parameters.getInputCost());
		graph_.insertEdge(v,t_);
		if(useZeroInOut) edgeScore.push_back(0);
		else edgeScore.push_back(parameters.getOutputCost());
	}


	for (int v = minVertex; v < maxVertex; ++v) {
		for (int i = 0; i < completeGraph.numberOfEdgesFromVertex(v); ++i) {
			size_t w=completeGraph.vertexFromVertex(v,i);
			if(w>maxVertex) continue;
			size_t e=completeGraph.edgeFromVertex(v,i);
			graph_.insertEdge(v-minVertex, w-minVertex);
			edgeScore.push_back(completeScore[e]);
		}
	}
	minV=minVertex;
	maxV=maxVertex;

}


//void LdpInstance::readGraph(std::ifstream& data,size_t maxVertex,char delim){
//	std::string line;
//	//	char delim = ' ';
//	size_t lineCounter=0;
//	std::getline(data, line);
//	lineCounter++;
//    if(debug()) std::cout << "called read graph" << std::endl;
//    //parameters.infoFile()<<"called read graph" << std::endl;
//	std::vector<std::string> strings = split(line, delim);
//	size_t numberOfVertices;

//	if (strings.size() == 1) {
//		if(parameters.isRestrictFrames()){
//			numberOfVertices=maxVertex+3;
//		}
//		else{
//			numberOfVertices = stoul(strings[0]);
//			numberOfVertices += 2;  //to include s and t in the end of the list
//		}
//		s_ = numberOfVertices - 2;
//		t_ = numberOfVertices - 1;

//	} else {
//		std::string str="first row must contain 1 number, detected ";
//		str+=std::to_string(strings.size());
//		str+="numbers";
//		throw std::runtime_error(str);
//	}



//	graphLifted_ = andres::graph::Digraph<>(numberOfVertices);
//	graph_ = andres::graph::Digraph<>(numberOfVertices);
//	std::vector<double> inputCosts(numberOfVertices-2,parameters.getInputCost());
//	std::vector<double> outputCosts(numberOfVertices-2,parameters.getOutputCost());


//	// std::vector<std::pair<size_t,siz  Data<>e_t> > liftedEdges;
//	vertexScore = std::vector<double>(numberOfVertices, 0);

//    if(debug())std::cout<<"Reading vertices from file. "<<std::endl;
////	parameters.infoFile()<<"Reading vertices from file. "<<std::endl;
////	parameters.infoFile().flush();
//	//Vertices that are not found have score=0. Appearance and disappearance cost are read here.
//	while (std::getline(data, line) && !line.empty()) {
//		lineCounter++;
//		strings = split(line, delim);
//		if (strings.size() < 2) {
//			throw std::runtime_error(
//					std::string("Vertex and its score expected"));
//		}


//		unsigned int v = std::stoul(strings[0]);
//		if(v>graph_.numberOfVertices()-3) continue;
//		double score = std::stod(strings[1]);
//		vertexScore[v] = score;

//		if(strings.size()==4){
//			inputCosts[v]=std::stod(strings[2]);
//			outputCosts[v]=std::stod(strings[3]);
//		}

//	}

//	for (int v = 0; v < numberOfVertices-2; ++v) {
//		graph_.insertEdge(s_,v);
//		edgeScore.push_back(inputCosts[v]);
//		graph_.insertEdge(v,t_);
//		edgeScore.push_back(outputCosts[v]);
//	}

//	size_t maxGap=parameters.getMaxTimeGapComplete();

//    if(debug()) std::cout<<"Reading base edges from file. "<<std::endl;
////	parameters.infoFile()<<"Reading base edges from file. "<<std::endl;
////	parameters.infoFile().flush();
//	while (std::getline(data, line) && !line.empty()) {
//		lineCounter++;
//		strings = split(line, delim);
//		if (strings.size() < 3) {
//			throw std::runtime_error(
//					std::string("Edge vertices and score expected, line "+std::to_string(lineCounter)));
//		}

//		unsigned int v = std::stoul(strings[0]);
//		unsigned int w = std::stoul(strings[1]);

//        if(v>numberOfVertices-3) break;
//        if(w>numberOfVertices-3) continue;

//		size_t gv=vertexGroups.getGroupIndex(v);
//		size_t gw=vertexGroups.getGroupIndex(w);
//		if(gw-gv>maxGap) continue;

//		//if(v>=graph_.numberOfVertices()-2||w>=graph_.numberOfVertices()-2) continue;
//		double score = std::stod(strings[2]);
//		auto edgeTest=graph_.findEdge(v,w);

//		if(!edgeTest.first){  //if the edge does not exist
//			graph_.insertEdge(v, w);
//			edgeScore.push_back(score);

//		}
//		else{  //if the edge already exists, only update the score
//			edgeScore[edgeTest.second]=score;

//		}
//	}
//    if(debug())std::cout<<"reading finished"<<std::endl;
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
        for (size_t j = 0; j < graphLifted_.numberOfEdgesFromVertex(i); ++j) {
           size_t vertex2=graphLifted_.vertexFromVertex(i,j);
           if(vertex2<minVertex){
               minVertex=vertex2;
           }
           if(vertex2>maxVertex){
               maxVertex=vertex2;
           }
        }
        liftedStructure[i]=ShiftedVector<char>(minVertex,maxVertex,false);
        for (size_t j = 0; j < graphLifted_.numberOfEdgesFromVertex(i); ++j) {
           size_t vertex2=graphLifted_.vertexFromVertex(i,j);
           liftedStructure[i][vertex2]=1;
        }
    }
}




void LdpInstance::sparsifyLiftedGraph(){


    parameters.getControlOutput()<<"Sparsify lifted graph"<<std::endl;
    parameters.writeControlOutput();
    //TODO run automaticLifted to find candidates first

    std::vector<double> newLiftedCosts;

    andres::graph::Digraph<> tempGraphLifted=(graph_.numberOfVertices());


    std::unordered_map<size_t,std::set<size_t>> liftedEdges;
    for (size_t v = 0; v < graphLifted_.numberOfVertices()-2; ++v) {
        std::unordered_set<size_t> alternativePath;
        for (size_t i = 0; i < graph_.numberOfEdgesFromVertex(v); ++i) {
            size_t w=graph_.vertexFromVertex(v,i);
            for(size_t u:reachable[w]){
                if(u!=w){
                    alternativePath.insert(u);
                }
            }
        }
        for (size_t i = 0; i < graph_.numberOfEdgesFromVertex(v); ++i) {
            size_t w=graph_.vertexFromVertex(v,i);
            if(alternativePath.count(w)==0){
                strongBaseEdges.at(v).insert(w);
            }
        }

        for (size_t i = 0; i < graphLifted_.numberOfEdgesFromVertex(v); ++i) {
            size_t w=graphLifted_.vertexFromVertex(v,i);
            if(w!=t_){
                if(alternativePath.count(w)>0) liftedEdges[v].insert(w);

            }
        }
    }



    parameters.getControlOutput()<<"done"<<std::endl;
   parameters.writeControlOutput();


    for (int i = 0; i < graphLifted_.numberOfEdges(); ++i) {
        size_t v0=graphLifted_.vertexOfEdge(i,0);
        size_t v1=graphLifted_.vertexOfEdge(i,1);
        int l0=vertexGroups.getGroupIndex(v0);
        int l1=vertexGroups.getGroupIndex(v1);
        double cost=getLiftedEdgeScore(i);
        bool goodCost=(cost<negativeLiftedThreshold)||(cost>positiveLiftedThreshold);
        if(isReachable(v0,v1)){

            int timeGapDiff=l1-l0-parameters.getDenseTimeLifted();
            bool timeConstraint=l1-l0<=parameters.getDenseTimeLifted()||((l1-l0)<=parameters.getMaxTimeLifted()&&(timeGapDiff%parameters.getLongerIntervalLifted())==0);
            if(timeConstraint&&goodCost){
                if(liftedEdges[v0].count(v1)>0){
                    tempGraphLifted.insertEdge(v0,v1);
                    newLiftedCosts.push_back(cost);
                }
                else{
                    auto edgeTest=graph_.findEdge(v0,v1);
                    if(edgeTest.first){
                        edgeScore[edgeTest.second]+=cost;  //Compensate that the lifted edge has been removed
                    }

                }

            }
        }

    }



    liftedEdgeScore=newLiftedCosts;

    graphLifted_=tempGraphLifted;
      parameters.getControlOutput()<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
    parameters.writeControlOutput();
    if(graphLifted_.numberOfEdges()!=newLiftedCosts.size()){
         parameters.getControlOutput()<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
        parameters.getControlOutput()<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
    }
    else{
         parameters.getControlOutput()<<"lifted edge number and lifted graph size match "<<std::endl;
        parameters.getControlOutput()<<"lifted edge number and lifted graph size match "<<std::endl;

    }
   parameters.writeControlOutput();
    initLiftedStructure();

}


//It could get complete structure and directly output two dim array
void LdpInstance::sparsifyBaseGraphNew(andres::graph::Digraph<> &inputGraph){
     parameters.getControlOutput()<<"Sparsify base graph"<<std::endl;
     parameters.writeControlOutput();

     std::vector<double> newBaseCosts;
     //std::vector<size_t> inOutEdges;
     size_t k=parameters.getKnnK();
     //std::vector<size_t> goodLongEdges;

   //  andres::graph::Digraph<> tempGraph(graph_.numberOfVertices());


     size_t nrVertices=inputGraph.numberOfVertices();
     std::vector<std::array<size_t,2>> edgesForMyGraph;

     std::vector<bool> finalEdges(inputGraph.numberOfEdges(),false);
     for (int v0 = 0; v0 < nrVertices; ++v0) {
         std::unordered_map<int,std::list<size_t>> edgesToKeep;
         size_t l0=vertexGroups.getGroupIndex(v0);
         for (size_t ne = 0; ne < inputGraph.numberOfEdgesFromVertex(v0); ++ne) {
             size_t e=inputGraph.edgeFromVertex(v0,ne);
             size_t v1=inputGraph.vertexFromVertex(v0,ne);
 //			std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<problemGraph.getEdgeScore(e)<<std::endl;
             if(v0==s_||v1==t_){
                 //tempGraph.insertEdge(v0,v1);
                 //newBaseCosts.push_back(edgeScore[e]);
                 finalEdges[e]=true;
             }
             else{
                 size_t l1=vertexGroups.getGroupIndex(v1);
                 size_t gap=l1-l0;
                 if(gap<=parameters.getMaxTimeBase()){
                 //if(gap<=parameters.getKnnTimeGap()){
                     //gap=std::min(parameters.getKnnTimeGap()+1,gap);
                     double cost=edgeScore[e];
                     if(edgesToKeep.count(gap)>0){
                         std::list<size_t>& smallList=edgesToKeep[gap];
                         auto it=smallList.begin();
                         double bsf=edgeScore[*it];
                         //std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
                         while(bsf>cost&&it!=smallList.end()){
                             it++;
                             size_t index=*it;
                             if(it!=smallList.end()){
                                 bsf=edgeScore[index];
                                 //	std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
                             }
                         }
                         if(it!=smallList.begin()){
                             smallList.insert(it,e);
                             if(smallList.size()>k) smallList.pop_front();
                         }
                         else if(smallList.size()<k){
                             smallList.push_front(e);
                         }
                     }
                     else{
                         edgesToKeep[gap].push_front(e);
                     }
                 }
 //				else if(gap<=parameters.getMaxTimeBase()){
 //					if(getEdgeScore(e)<=parameters.getBaseUpperThreshold()){
 //						//tempGraph.insertEdge(v0,v1);
 //						//newBaseCosts.push_back(getEdgeScore(e));
 //						finalEdges[e]=true;
 //					}
 //
 //				}
             }
         }
         //std::cout.precision(4);
         double bsf=0;
         for (int gap = 0; gap <= parameters.getKnnTimeGap(); ++gap) {
             if(edgesToKeep.count(gap)>0){
                 auto& smallList=edgesToKeep[gap];
                 for(size_t e:smallList){
                     finalEdges[e]=true;
                     if(edgeScore[e]<bsf){
                         bsf=edgeScore[e];
                     }
                 }
             }
         }

         for (int gap =  parameters.getKnnTimeGap()+1;gap<=parameters.getMaxTimeBase(); ++gap) {
             if(edgesToKeep.count(gap)>0){
                 auto& smallList=edgesToKeep[gap];
                 for(size_t e:smallList){
                     double score=edgeScore[e];
                     if(score<=baseThreshold){
                         finalEdges[e]=true;
                     }
                 }

             }
         }

     }

     for (int e = 0; e < inputGraph.numberOfEdges(); ++e) {
         if(finalEdges[e]){
             size_t v0=inputGraph.vertexOfEdge(e,0);
             size_t v1=inputGraph.vertexOfEdge(e,1);
             //tempGraph.insertEdge(v0,v1);
             edgesForMyGraph.push_back({v0,v1});
             newBaseCosts.push_back(edgeScore[e]);
         }
     }

     assert(edgesForMyGraph.size()==newBaseCosts.size());

     edgeScore=newBaseCosts;


     if(edgesForMyGraph.size()!=newBaseCosts.size()){
         parameters.getControlOutput()<<"edge number mismatch, edge vector: "<<edgesForMyGraph.size()<<", cost vector "<<newBaseCosts.size()<<std::endl;
         parameters.writeControlOutput();

     }
     else{

         parameters.getControlOutput()<<"edge number and graph size match "<<std::endl;
         parameters.writeControlOutput();
     }

     parameters.getControlOutput()<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
     parameters.writeControlOutput();

     myGraph=LdpDirectedGraph(edgesForMyGraph,newBaseCosts);

 }



//It could get complete structure and directly output two dim array
void LdpInstance::sparsifyBaseGraph(){
     parameters.getControlOutput()<<"Sparsify base graph"<<std::endl;
     parameters.writeControlOutput();

     std::vector<double> newBaseCosts;
     //std::vector<size_t> inOutEdges;
     size_t k=parameters.getKnnK();
     //std::vector<size_t> goodLongEdges;

     andres::graph::Digraph<> tempGraph(graph_.numberOfVertices());



     std::vector<bool> finalEdges(graph_.numberOfEdges(),false);
     for (int v0 = 0; v0 < graph_.numberOfVertices(); ++v0) {
         std::unordered_map<int,std::list<size_t>> edgesToKeep;
         size_t l0=vertexGroups.getGroupIndex(v0);
         for (size_t ne = 0; ne < graph_.numberOfEdgesFromVertex(v0); ++ne) {
             size_t e=graph_.edgeFromVertex(v0,ne);
             size_t v1=graph_.vertexFromVertex(v0,ne);
 //			std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<problemGraph.getEdgeScore(e)<<std::endl;
             if(v0==s_||v1==t_){
                 //tempGraph.insertEdge(v0,v1);
                 //newBaseCosts.push_back(edgeScore[e]);
                 finalEdges[e]=true;
             }
             else{
                 size_t l1=vertexGroups.getGroupIndex(v1);
                 size_t gap=l1-l0;
                 if(gap<=parameters.getMaxTimeBase()){
                 //if(gap<=parameters.getKnnTimeGap()){
                     //gap=std::min(parameters.getKnnTimeGap()+1,gap);
                     double cost=edgeScore[e];
                     if(edgesToKeep.count(gap)>0){
                         std::list<size_t>& smallList=edgesToKeep[gap];
                         auto it=smallList.begin();
                         double bsf=edgeScore[*it];
                         //std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
                         while(bsf>cost&&it!=smallList.end()){
                             it++;
                             size_t index=*it;
                             if(it!=smallList.end()){
                                 bsf=edgeScore[index];
                                 //	std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<bsf<<std::endl;
                             }
                         }
                         if(it!=smallList.begin()){
                             smallList.insert(it,e);
                             if(smallList.size()>k) smallList.pop_front();
                         }
                         else if(smallList.size()<k){
                             smallList.push_front(e);
                         }
                     }
                     else{
                         edgesToKeep[gap].push_front(e);
                     }
                 }
 //				else if(gap<=parameters.getMaxTimeBase()){
 //					if(getEdgeScore(e)<=parameters.getBaseUpperThreshold()){
 //						//tempGraph.insertEdge(v0,v1);
 //						//newBaseCosts.push_back(getEdgeScore(e));
 //						finalEdges[e]=true;
 //					}
 //
 //				}
             }
         }
         //std::cout.precision(4);
         double bsf=0;
         for (int gap = 0; gap <= parameters.getKnnTimeGap(); ++gap) {
             if(edgesToKeep.count(gap)>0){
                 auto& smallList=edgesToKeep[gap];
                 for(size_t e:smallList){
                     finalEdges[e]=true;
                     if(edgeScore[e]<bsf){
                         bsf=edgeScore[e];
                     }
                 }
             }
         }

         for (int gap =  parameters.getKnnTimeGap()+1;gap<=parameters.getMaxTimeBase(); ++gap) {
             if(edgesToKeep.count(gap)>0){
                 auto& smallList=edgesToKeep[gap];
                 for(size_t e:smallList){
                     double score=edgeScore[e];
                     if(score<=baseThreshold){
                         finalEdges[e]=true;
                     }
                 }

             }
         }

     }

     for (int e = 0; e < graph_.numberOfEdges(); ++e) {
         if(finalEdges[e]){
             size_t v0=graph_.vertexOfEdge(e,0);
             size_t v1=graph_.vertexOfEdge(e,1);
             tempGraph.insertEdge(v0,v1);
             newBaseCosts.push_back(edgeScore[e]);
         }
     }

     if(newBaseCosts.size()!=tempGraph.numberOfEdges()){
         throw std::runtime_error("Error in base graph sparsification.");
     }



     graph_=tempGraph;
     edgeScore=newBaseCosts;


     if(graph_.numberOfEdges()!=newBaseCosts.size()){
         parameters.getControlOutput()<<"edge number mismatch, graph: "<<graph_.numberOfEdges()<<", cost vector "<<newBaseCosts.size()<<std::endl;
         parameters.writeControlOutput();

     }
     else{

         parameters.getControlOutput()<<"edge number and graph size match "<<std::endl;
         parameters.writeControlOutput();
     }

     parameters.getControlOutput()<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
     parameters.writeControlOutput();

 }










bool LdpInstance::isStrongBase(size_t v,size_t w) const{
    bool isStrong= strongBaseEdges.at(v).count(w)>0;
    return isStrong;
   // return false;
}








}}//End of namespaces
