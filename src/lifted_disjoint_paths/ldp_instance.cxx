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
//    disjointPaths::CompleteStructure<> csBase(*twoGraphsIS.myPvg);
//    //TODO pointer to VG cannot be constant!
//    csBase.addEdgesFromVectorsAll(*twoGraphsIS.pBaseEdges,*twoGraphsIS.pBaseCosts);
//    graph_=csBase.completeGraph;
//    edgeScore=csBase.completeScore;
//    //TODO base graph needs s,t edges and scores!

//    disjointPaths::CompleteStructure<> csLifted(*twoGraphsIS.myPvg);
//    csLifted.addEdgesFromVectorsAll(*twoGraphsIS.pLiftedEdges,*twoGraphsIS.pLiftedCosts);
//    graphLifted_=csLifted.completeGraph;
//    liftedEdgeScore=csLifted.completeScore;

//    //TODO set all remaining global variables
//}

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
        disjointPaths::keepFractionOfLifted(*this,parameters);
        initLiftedStructure();
      //  sparsifyLiftedGraph();

    }
    else{
        //desc=initReachable(graph_,parameters);
        reachable=disjointPaths::initReachableSet(graph_,parameters,&vertexGroups);
        initLiftedStructure();
    }

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
    liftedStructure=std::vector<ShiftedVector<bool>>(graphLifted_.numberOfVertices());
    for (size_t i = 0; i < graphLifted_.numberOfVertices(); ++i) {
        size_t maxVertex=0;
        size_t minVertex=graphLifted_.numberOfVertices();
        for (size_t j = 0; j < graphLifted_.numberOfEdgesFromVertex(i); ++j) {
           size_t vertex2=graphLifted_.vertexFromVertex(i,j);
           if(vertex2<minVertex){
               minVertex=vertex2;
           }
           if(vertex2>maxVertex){
               maxVertex=vertex2;
           }
        }
        liftedStructure[i]=ShiftedVector<bool>(minVertex,maxVertex,false);
        for (size_t j = 0; j < graphLifted_.numberOfEdgesFromVertex(i); ++j) {
           size_t vertex2=graphLifted_.vertexFromVertex(i,j);
           liftedStructure[i].setValue(vertex2,true);
        }
    }
}





bool LdpInstance::isStrongBase(size_t v,size_t w) const{
    bool isStrong= strongBaseEdges.at(v).count(w)>0;
    return isStrong;
   // return false;
}



}}//End of namespaces
