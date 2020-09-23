#include"lifted_disjoint_paths/ldp_instance.hxx"


namespace LPMP{
namespace lifted_disjoint_paths {






LdpInstance::LdpInstance(const ConfigDisjoint<> &configParameters, disjointPaths::CompleteStructure<>& cs):
    parameters(configParameters)
{

    size_t minT=1;
    size_t maxT=cs.maxTime+1;
    readGraphWithTime(minT,maxT,&cs);

    init();

}








LdpInstance::LdpInstance(const ConfigDisjoint<>& configParameters):
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

void LdpInstance::init(){

    std::cout<<"Adding automatic lifted edges"<<std::endl;
    parameters.infoFile()<<"Adding automatic lifted edges"<<std::endl;
    parameters.infoFile().flush();
    for (size_t i = 0; i < graph_.numberOfEdges(); ++i) {
        size_t v0=graph_.vertexOfEdge(i,0);
        size_t v1=graph_.vertexOfEdge(i,1);
        if(v0!=s_&&v1!=t_){
            //	if(secOrderDesc[v0][v1]){
            graphLifted_.insertEdge(v0,v1);
            liftedEdgeScore.push_back(edgeScore[i]);

        }
    }
    std::cout<<"done"<<std::endl;
    parameters.infoFile()<<"done"<<std::endl;
    parameters.infoFile().flush();


    std::cout<<"number of vertices "<<graph_.numberOfVertices()<<std::endl;
    parameters.infoFile()<<"number of vertices "<<graph_.numberOfVertices()<<std::endl;
    parameters.infoFile().flush();
    strongBaseEdges=std::vector<std::unordered_set<size_t>>(graph_.numberOfVertices());
    if(parameters.isSparsify()){
        sparsifyBaseGraph();
        sparsifyLiftedGraph();

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
	std::cout << "called read graph" << std::endl;
	parameters.infoFile()<<"called read graph" << std::endl;
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

	std::cout<<"Reading vertices from file. "<<std::endl;
	parameters.infoFile()<<"Reading vertices from file. "<<std::endl;
	parameters.infoFile().flush();
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

	std::cout<<"Reading base edges from file. "<<std::endl;
	parameters.infoFile()<<"Reading base edges from file. "<<std::endl;
	parameters.infoFile().flush();
	while (std::getline(data, line) && !line.empty()) {
		lineCounter++;
		strings = split(line, delim);
		if (strings.size() < 3) {
			throw std::runtime_error(
					std::string("Edge vertices and score expected, line "+std::to_string(lineCounter)));
		}

		unsigned int v = std::stoul(strings[0]);
		unsigned int w = std::stoul(strings[1]);

		if(v>numberOfVertices-3||w>numberOfVertices-3) continue;

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
    std::cout<<"reading finished"<<std::endl;
}


//template<typename LABEL_ITERATOR>
//bool LdpInstance::check_feasiblity(LABEL_ITERATOR begin, LABEL_ITERATOR end) const{


//	//Labels into vector
//	size_t indexCounter=0;
//	std::vector<double> solution0(graph_.numberOfEdges()+graphLifted_.numberOfEdges()+graph_.numberOfVertices());
//	auto it=begin;
//	for (;it!=end&&indexCounter<solution0.size();it++) {
//		double value=*it;
//		solution0[indexCounter]=value;
//		indexCounter++;
//	}

//	if(it!=end){
//		throw std::runtime_error("Wrong size of solution labels for feasibility check.");
//	}


//	//Check flow conservation
//	bool isFeasible=true;
//	std::unordered_map<size_t,size_t> predecessors; //Used later for lifted edges consistency
//	for (int i = 0; i < graph_.numberOfVertices()&&isFeasible; ++i) {
//		if(i==s_||i==t_) continue;
//		bool vertexActive=solution0[getVertexVarIndex(i)]>0.5;
//		bool outputFound=false;
//		for (int j = 0; j < graph_.numberOfEdgesFromVertex(i)&&isFeasible; ++j) {
//			size_t e=graph_.edgeFromVertex(i,j);
//			if(solution0[getEdgeVarIndex(e)]>0.5){
//				if(outputFound||!vertexActive){
//					isFeasible=false;
//				}

//				else{
//					outputFound=true;
//				}
//			}
//		}
//		if(vertexActive&&!outputFound){
//			isFeasible=false;
//		}

//		bool inputFound=false;
//		for (int j = 0; j < graph_.numberOfEdgesToVertex(i)&&isFeasible; ++j) {
//			size_t e=graph_.edgeToVertex(i,j);
//			if(solution0[getEdgeVarIndex(e)]>0.5){
//				if(inputFound||!vertexActive){
//					isFeasible=false;
//				}
//				else{
//					inputFound=true;
//					predecessors[i]=graph_.vertexToVertex(i,j);
//				}
//			}
//		}
//		if(vertexActive&&!inputFound){
//			isFeasible=false;
//		}

//	}

//	//Check lifted edge labels consistency with paths
//	for (int i = 0; i < graph_.numberOfEdgesToVertex(t_)&&isFeasible; ++i) {
//		size_t e=graph_.edgeToVertex(t_,i);
//		size_t eVarIndex=getEdgeVarIndex(e);
//		if(solution0[eVarIndex] > 0.5){

//			std::vector<bool> isOnPath(graph_.numberOfVertices(),0);
//			size_t vertex=graph_.vertexToVertex(t_,i);
//			isOnPath[vertex]=1;
//			isOnPath[t_]=1;  //Maybe not necessary

//			while(vertex!=s_&&isFeasible){

//				for (int j = 0; j < graphLifted_.numberOfEdgesFromVertex(vertex); ++j) {
//					size_t le=graphLifted_.edgeFromVertex(vertex,j);
//					size_t vertex2=graphLifted_.vertexFromVertex(vertex,j);
//					size_t leVarIndex=getLiftedEdgeVarIndex(le);
//					//size_t vertex2VarIndex=data_.getVertexVarIndex(vertex2);
//					if(isOnPath[vertex2]&&solution0[leVarIndex]<0.5) {
//						isFeasible=false;
//					}
//					else if(!isOnPath[vertex2]&&solution0[leVarIndex]>0.5){
//						isFeasible=false;
//					}
//				}
//				isOnPath[vertex]=1;
//				vertex=predecessors[vertex];
//			}

//		}
//	}

//	return isFeasible;

//}



//template<typename EDGE_LABEL_ITERATOR>
//double LdpInstance::evaluate(EDGE_LABEL_ITERATOR begin, EDGE_LABEL_ITERATOR end) const{
//	size_t indexCounter=0;
//	std::vector<double> solution0(graph_.numberOfEdges()+graphLifted_.numberOfEdges()+graph_.numberOfVertices());
//	auto it=begin;
//	for (;it!=end&&indexCounter<solution0.size();it++) {
//		double value=*it;
//		solution0[indexCounter]=value;
//		indexCounter++;
//	}

//	if(it!=end){
//		throw std::runtime_error("Wrong size of solution labels for evaluation.");
//	}

//	double objectiveValue=0;
//	for (int i = 0; i < numberOfEdges; ++i) {
//		if(solution0[getEdgeVarIndex(i)]>0.5) objectiveValue+=edgeScore[i];
//	}
//	for (int i = 0; i < numberOfLiftedEdges; ++i) {
//		if(solution0[getLiftedEdgeVarIndex(i)]>0.5) objectiveValue+=liftedEdgeScore[i];
//	}

//	return objectiveValue;


//}

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



void LdpInstance::sparsifyBaseGraph(){
	std::cout<<"Sparsify base graph"<<std::endl;
	parameters.infoFile()<<"Sparsify base graph"<<std::endl;
	andres::graph::Digraph<> tempGraph(graph_.numberOfVertices());
	std::vector<double> newBaseCosts;
	//std::vector<size_t> inOutEdges;
	size_t k=parameters.getKnnK();
	//std::vector<size_t> goodLongEdges;


	std::vector<bool> finalEdges(graph_.numberOfEdges(),false);
	for (int v0 = 0; v0 < graph_.numberOfVertices(); ++v0) {
		std::unordered_map<int,std::list<size_t>> edgesToKeep;
		size_t l0=vertexGroups.getGroupIndex(v0);
		for (size_t ne = 0; ne < graph_.numberOfEdgesFromVertex(v0); ++ne) {
			size_t e=graph_.edgeFromVertex(v0,ne);
			size_t v1=graph_.vertexFromVertex(v0,ne);
			//			std::cout<<"edge "<<e<<": "<<v0<<","<<v1<<": "<<problemGraph.getEdgeScore(e)<<std::endl;
			if(v0==s_||v1==t_){

				finalEdges[e]=true;
			}
			else{
				size_t l1=vertexGroups.getGroupIndex(v1);
				size_t gap=l1-l0;
				if(gap<=parameters.getMaxTimeBase()){

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

			}
		}

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
					if(score<=parameters.getBaseUpperThreshold()){
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
		parameters.infoFile()<<"edge number mismatch, graph: "<<graph_.numberOfEdges()<<", cost vector "<<newBaseCosts.size()<<std::endl;
		parameters.infoFile().flush();
		std::cout<<"edge number mismatch, graph: "<<graph_.numberOfEdges()<<", cost vector "<<newBaseCosts.size()<<std::endl;
	}
	else{
		std::cout<<"edge number and graph size match "<<std::endl;
		parameters.infoFile()<<"edge number and graph size match "<<std::endl;
		parameters.infoFile().flush();
	}

    reachable=initReachableSet(graph_,parameters,&vertexGroups);


	std::cout<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
	parameters.infoFile()<<"Left "<<newBaseCosts.size()<<" base edges"<<std::endl;
	parameters.infoFile().flush();

}

void LdpInstance::sparsifyLiftedGraph(){

	std::cout<<"Sparsify lifted graph"<<std::endl;
	parameters.infoFile()<<"Sparsify lifted graph"<<std::endl;
	parameters.infoFile().flush();
	//TODO run automaticLifted to find candidates first

	double negMaxValue=0;
	double posMinValue=0;
	bool useAdaptive=false;
	//TODO adaptive lifted threshold in config file
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


	std::cout<<"done"<<std::endl;
	parameters.infoFile()<<"done"<<std::endl;
	parameters.infoFile().flush();


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
	std::cout<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
	parameters.infoFile()<<"Left "<<newLiftedCosts.size()<<" lifted edges."<<std::endl;
	parameters.infoFile().flush();

	if(graphLifted_.numberOfEdges()!=newLiftedCosts.size()){
		std::cout<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
		parameters.infoFile()<<"lifted edge number mismatch, lifted graph: "<<graphLifted_.numberOfEdges()<<", cost vector "<<newLiftedCosts.size()<<std::endl;
	}
	else{
		std::cout<<"lifted edge number and lifted graph size match "<<std::endl;
		parameters.infoFile()<<"lifted edge number and lifted graph size match "<<std::endl;

	}
	parameters.infoFile().flush();
    initLiftedStructure();

}

bool LdpInstance::isStrongBase(size_t v,size_t w) const{
    bool isStrong= strongBaseEdges.at(v).count(w)>0;
    return isStrong;
   // return false;
}



}}//End of namespaces
