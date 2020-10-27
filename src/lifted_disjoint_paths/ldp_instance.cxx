#include"lifted_disjoint_paths/ldp_instance.hxx"


namespace LPMP{
namespace lifted_disjoint_paths {






LdpInstance::LdpInstance(LdpParameters<> &configParameters, disjointPaths::CompleteStructure<>& cs):
    parameters(configParameters)
{

    size_t minT=1;
    size_t maxT=cs.maxTime+1;
    readGraphWithTime(minT,maxT,&cs);

    init();

}








LdpInstance::LdpInstance( LdpParameters<>& configParameters):
    parameters(configParameters)
{
    char delim=',';
    size_t maxVertex;

        vertexGroups=disjointPaths::VertexGroups<size_t>(parameters);
        maxVertex=vertexGroups.getMaxVertex();

    std::ifstream graphFile(parameters.getGraphFileName());
    readGraph(graphFile,maxVertex,delim);
    graphFile.close();

    init();

    //baseEdgeLabels=std::vector<bool>(numberOfEdges);
}



//LdpInstance::LdpInstance(LdpParameters<>& configParameters,const disjointPaths::TwoGraphsInputStructure& twoGraphsIS):parameters(configParameters){
LdpInstance::LdpInstance(LdpParameters<>& configParameters,const py::array_t<size_t>& baseEdges,const py::array_t<size_t>& liftedEdges,const  py::array_t<double>& baseCosts,const  py::array_t<double>& liftedCosts,disjointPaths::VertexGroups<>& pvg):parameters(configParameters){
//LdpInstance::LdpInstance(LdpParameters<>& configParameters,const std::vector<std::array<size_t,2>>& baseEdges,const std::vector<std::array<size_t,2>>& liftedEdges,const  std::vector<double>& baseCosts,const  std::vector<double>& liftedCosts,disjointPaths::VertexGroups<>& pvg):parameters(configParameters){
    disjointPaths::CompleteStructure<> csBase(pvg);

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

    disjointPaths::CompleteStructure<> csLifted(vertexGroups);
    csLifted.addEdgesFromVectorsAll(liftedEdges,liftedCosts);
    graphLifted_=csLifted.completeGraph;
    liftedEdgeScore=csLifted.completeScore;


    if(configParameters.isSparsify()){
        disjointPaths::createKnnBaseGraph(*this,configParameters);
        reachable=disjointPaths::initReachableSet(graph_,parameters,&vertexGroups);
        //disjointPaths::keepFractionOfLifted(*this,configParameters);
        sparsifyLiftedGraph();
    }
    else{
        reachable=disjointPaths::initReachableSet(graph_,parameters,&vertexGroups);
    }



    initLiftedStructure();

    numberOfVertices=graph_.numberOfVertices();
    numberOfLiftedEdges=graphLifted_.numberOfEdges();
    numberOfEdges=graph_.numberOfEdges();

    vertexScore=std::vector<double>(graph_.numberOfVertices()-2);

    myGraph=LdpDirectedGraph(graph_,edgeScore);
    myGraphLifted=LdpDirectedGraph(graphLifted_,liftedEdgeScore);

    minV=0;
    maxV=maxVertex;


}

void LdpInstance::init(){

    if(debug()) std::cout<<"Adding automatic lifted edges"<<std::endl;
  //  parameters.infoFile()<<"Adding automatic lifted edges"<<std::endl;
   // parameters.infoFile().flush();
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
    if(debug()) std::cout<<"done"<<std::endl;
   // parameters.infoFile()<<"done"<<std::endl;
  //  parameters.infoFile().flush();


   if(diagnostics())  std::cout<<"number of vertices "<<graph_.numberOfVertices()<<std::endl;
  //  parameters.infoFile()<<"number of vertices "<<graph_.numberOfVertices()<<std::endl;
  //  parameters.infoFile().flush();
    strongBaseEdges=std::vector<std::unordered_set<size_t>>(graph_.numberOfVertices());
    if(parameters.isSparsify()){
        //sparsifyBaseGraph();
        disjointPaths::createKnnBaseGraph(*this,parameters);
        reachable=initReachableSet(graph_,parameters,&vertexGroups);
       // disjointPaths::keepFractionOfLifted(*this,parameters);
        sparsifyLiftedGraph();
        initLiftedStructure();
      //  sparsifyLiftedGraph();well

    }
    else{
        //desc=initReachable(graph_,parameters);
        reachable=disjointPaths::initReachableSet(graph_,parameters,&vertexGroups);
        initLiftedStructure();
    }

    myGraph=LdpDirectedGraph(graph_,edgeScore);
    myGraphLifted=LdpDirectedGraph(graphLifted_,liftedEdgeScore);

    numberOfVertices=graph_.numberOfVertices();
    numberOfEdges=graph_.numberOfEdges();
    numberOfLiftedEdges=graphLifted_.numberOfEdges();

}



void LdpInstance::readGraphWithTime(size_t minTime,size_t maxTime,disjointPaths::CompleteStructure<>* cs){

	andres::graph::Digraph<>& completeGraph=cs->completeGraph;
	std::vector<double>& completeScore=cs->completeScore;
    disjointPaths::VertexGroups<> vg=cs->getVertexGroups();

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

	for (int gi = minTime; gi < maxTime; ++gi) {
		//groups[gi-minTime+1]=std::vector<size_t>();
		for(size_t v:vg.getGroupVertices(gi)){
			size_t vertex=v-minVertex;
			groups[gi-minTime+1].push_back(vertex);
			vToGroup[vertex]=gi-minTime+1;
		}
	}

    vertexGroups=disjointPaths::VertexGroups<>(groups,vToGroup);
	vertexScore = std::vector<double>(numberOfVertices, 0);
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


void LdpInstance::readGraph(std::ifstream& data,size_t maxVertex,char delim){
	std::string line;
	//	char delim = ' ';
	size_t lineCounter=0;
	std::getline(data, line);
	lineCounter++;
    if(debug()) std::cout << "called read graph" << std::endl;
    //parameters.infoFile()<<"called read graph" << std::endl;
	std::vector<std::string> strings = split(line, delim);
	size_t numberOfVertices;

	if (strings.size() == 1) {
		if(parameters.isRestrictFrames()){
			numberOfVertices=maxVertex+3;
		}
		else{
			numberOfVertices = stoul(strings[0]);
			numberOfVertices += 2;  //to include s and t in the end of the list
		}
		s_ = numberOfVertices - 2;
		t_ = numberOfVertices - 1;

	} else {
		std::string str="first row must contain 1 number, detected ";
		str+=std::to_string(strings.size());
		str+="numbers";
		throw std::runtime_error(str);
	}



	graphLifted_ = andres::graph::Digraph<>(numberOfVertices);
	graph_ = andres::graph::Digraph<>(numberOfVertices);
	std::vector<double> inputCosts(numberOfVertices-2,parameters.getInputCost());
	std::vector<double> outputCosts(numberOfVertices-2,parameters.getOutputCost());


	// std::vector<std::pair<size_t,siz  Data<>e_t> > liftedEdges;
	vertexScore = std::vector<double>(numberOfVertices, 0);

    if(debug())std::cout<<"Reading vertices from file. "<<std::endl;
//	parameters.infoFile()<<"Reading vertices from file. "<<std::endl;
//	parameters.infoFile().flush();
	//Vertices that are not found have score=0. Appearance and disappearance cost are read here.
	while (std::getline(data, line) && !line.empty()) {
		lineCounter++;
		strings = split(line, delim);
		if (strings.size() < 2) {
			throw std::runtime_error(
					std::string("Vertex and its score expected"));
		}


		unsigned int v = std::stoul(strings[0]);
		if(v>graph_.numberOfVertices()-3) continue;
		double score = std::stod(strings[1]);
		vertexScore[v] = score;

		if(strings.size()==4){
			inputCosts[v]=std::stod(strings[2]);
			outputCosts[v]=std::stod(strings[3]);
		}

	}

	for (int v = 0; v < numberOfVertices-2; ++v) {
		graph_.insertEdge(s_,v);
		edgeScore.push_back(inputCosts[v]);
		graph_.insertEdge(v,t_);
		edgeScore.push_back(outputCosts[v]);
	}

	size_t maxGap=parameters.getMaxTimeGapComplete();

    if(debug()) std::cout<<"Reading base edges from file. "<<std::endl;
//	parameters.infoFile()<<"Reading base edges from file. "<<std::endl;
//	parameters.infoFile().flush();
	while (std::getline(data, line) && !line.empty()) {
		lineCounter++;
		strings = split(line, delim);
		if (strings.size() < 3) {
			throw std::runtime_error(
					std::string("Edge vertices and score expected, line "+std::to_string(lineCounter)));
		}

		unsigned int v = std::stoul(strings[0]);
		unsigned int w = std::stoul(strings[1]);

        if(v>numberOfVertices-3) break;
        if(w>numberOfVertices-3) continue;

		size_t gv=vertexGroups.getGroupIndex(v);
		size_t gw=vertexGroups.getGroupIndex(w);
		if(gw-gv>maxGap) continue;

		//if(v>=graph_.numberOfVertices()-2||w>=graph_.numberOfVertices()-2) continue;
		double score = std::stod(strings[2]);
		auto edgeTest=graph_.findEdge(v,w);

		if(!edgeTest.first){  //if the edge does not exist
			graph_.insertEdge(v, w);
			edgeScore.push_back(score);

		}
		else{  //if the edge already exists, only update the score
			edgeScore[edgeTest.second]=score;

		}
	}
    if(debug())std::cout<<"reading finished"<<std::endl;
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

    double negMaxValue=0;
    double posMinValue=0;

    std::vector<double> newLiftedCosts;

    andres::graph::Digraph<> tempGraphLifted=(graph_.numberOfVertices());


    negMaxValue=parameters.getNegativeThresholdLifted();
    posMinValue=parameters.getPositiveThresholdLifted();


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
        bool goodCost=(cost<negMaxValue)||(cost>posMinValue);
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

bool LdpInstance::isStrongBase(size_t v,size_t w) const{
    bool isStrong= strongBaseEdges.at(v).count(w)>0;
    return isStrong;
   // return false;
}








}}//End of namespaces
