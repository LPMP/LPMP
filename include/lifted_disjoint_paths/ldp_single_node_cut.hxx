#pragma once
#include <andres/graph/digraph.hxx>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <array>
#include <list>


namespace LPMP {

//struct baseAndLiftedMessages{
//	std::unordered_map<size_t,double> baseMessages;
//	std::unordered_map<size_t,double> liftedMessages;
//};

struct StrForUpdateValues{
	std::unordered_map<size_t,double>& solutionCosts;
	mutable std::unordered_map<size_t,double> valuesStructure;
	mutable std::unordered_map<size_t,size_t> indexStructure;
	const std::unordered_map<size_t,double>& baseCosts;
	const std::unordered_map<size_t,double>& liftedCosts;

	//std::unordered_set<size_t> relevantVertices;
	bool useAllVertices;
	double optValue;
	const size_t nodeID;
	//size_t optimalSolution;
	StrForUpdateValues(const std::unordered_map<size_t,double>& bCosts,const std::unordered_map<size_t,double>& lCosts,std::unordered_map<size_t,double>& solCosts,const size_t centralNode):
	baseCosts(bCosts),
	liftedCosts(lCosts),
	solutionCosts(solCosts),
	optValue(0),
	nodeID(centralNode)
	{
		useAllVertices=true;
		//solutionCosts=baseCosts;
		//optimalSolution=nodeID;
	}
	bool useVertex(size_t vertex){
		if(useAllVertices) return true;
		//return relevantVertices.count(vertex)>0;
		return valuesStructure.count(vertex)>0;
	}

	void setUseAllVertices(bool value){
		useAllVertices=value;
	}
	double getOptimalValue(){
		return optValue;
	}
//	double getOptimalValue(){
//		auto it=indexStructure.find(nodeID);
//		if(it!=indexStructure.end()){
//			return solutionCosts.at(it->second);
//		}
//		else{
//			return 0;
//		}
//	}
};

template<class LDP_INSTANCE>
class ldp_single_node_cut_factor
{
public:
	//constexpr static std::size_t no_edge_active = std::numeric_limits<std::size_t>::infinity();

	//template<class LPD_STRUCT> ldp_single_node_cut_factor(const LPD_STRUCT& ldpStruct);
	ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut);

	void initBaseCosts(double fractionBase);

	void initLiftedCosts(double fractionLifted);

	void addLiftedEdge(size_t node,double cost);

	const andres::graph::Digraph<>& getBaseGraph() const {
		return baseGraph;
	}

	const std::unordered_map<size_t, double>& getLiftedCosts() const {
		return liftedCosts;
	}

	const std::unordered_map<size_t, double>& getBaseCosts() const {
		return baseCosts;
	}


	bool hasLiftedEdgeToNode(size_t vertex){
		return liftedCosts.count(vertex)>0;
	}

	double LowerBound() const;


	double EvaluatePrimal() const;

	void init_primal(){
		primalBase_=nodeNotActive;
	}

	void setBaseEdgeActive(size_t vertex);
	void setNoBaseEdgeActive(size_t vertex);

	void setPrimalLifted(std::unordered_set<size_t>& verticesOfActiveEdges);

	bool isNodeActive(){ return primalBase_!=nodeNotActive;}

	size_t getPrimalBase() const {
		return primalBase_;
	}

	const std::unordered_set<size_t>& getPrimalLifted() const {
		return primalLifted_;
	}


	template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
	template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

	//auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }//TODO change this. This will not work with so many variables
	double tmp_to_delete_val;
	auto export_variables() { return std::tie(tmp_to_delete_val); } //?What comes here?

	void updateCostSimple(const double value,const size_t vertexIndex,bool isLifted);
	double getOneBaseEdgeMinMarginal(size_t vertex);
	std::unordered_map<size_t,double> getAllBaseMinMarginals();

	std::unordered_map<size_t,double> getAllLiftedMinMarginals();

	double oneLiftedMinMarginal(size_t vertexOfLiftedEdge);



	const std::size_t nodeID;
	const size_t nodeNotActive;
	size_t primalBase_;
	std::unordered_set<size_t> primalLifted_;




private:
	void updateValues() const;
	void updateValues(StrForUpdateValues& myStr,size_t vertexToIgnore=std::numeric_limits<size_t>::max()) const;
	std::unordered_map<size_t,double> bottomUpUpdate(const StrForUpdateValues& myStr,size_t vertex,std::unordered_set<size_t>* pClosedVert=0,std::unordered_map<size_t,double>* pBUValuesStr=0);
	void updateOptimal() const;

	std::unordered_map<size_t,std::unordered_set<size_t>> createCandidateGraph(const StrForUpdateValues& myStr);


	void findAllOptimal(std::unordered_set<size_t>& isNotZeroInOpt,std::unordered_set<size_t>& isOneInOpt,std::unordered_map<size_t,std::unordered_set<size_t>>& candidateGraph,const StrForUpdateValues& strForUpdateValues);

	size_t getNeighborBaseVertex(size_t firstNode,size_t neighborIndex) const{
		assert(firstNode < baseGraph.numberOfVertices());
		if(isOutFlow){
			return baseGraph.vertexFromVertex(firstNode,neighborIndex);
		}
		else{
			return baseGraph.vertexToVertex(firstNode,neighborIndex);
		}
	}
	size_t numberOfNeighborsBase(const size_t nodeIndex) const {
		assert(nodeIndex < baseGraph.numberOfVertices());
		if(isOutFlow){
			return baseGraph.numberOfEdgesFromVertex(nodeIndex);
		}
		else{
			return baseGraph.numberOfEdgesToVertex(nodeIndex);
		}
	}
	size_t numberOfNeighborsBaseRev(const size_t nodeIndex) const {
		assert(nodeIndex < baseGraph.numberOfVertices());
		if(!isOutFlow){
			return baseGraph.numberOfEdgesFromVertex(nodeIndex);
		}
		else{
			return baseGraph.numberOfEdgesToVertex(nodeIndex);
		}
	}
	bool isInThisFactorRange(const size_t nodeIndex) const {
		assert(nodeIndex < baseGraph.numberOfVertices());
		if(isOutFlow){
			if(nodeIndex==ldpInstance.getTerminalNode()) return true;
			else return ldpInstance.getGroupIndex(nodeIndex)<=maxLayer;
		}
		else{
			if(nodeIndex==ldpInstance.getSourceNode()) return true;
			else return ldpInstance.getGroupIndex(nodeIndex)>=minLayer;
		}
	}

	bool isInGivenRange(const size_t nodeIndex,const size_t boundLayer) const {
		assert(nodeIndex < baseGraph.numberOfVertices());
		if(isOutFlow){
			return ldpInstance.getGroupIndex(nodeIndex)<=boundLayer;
		}
		else{
			return ldpInstance.getGroupIndex(nodeIndex)>=boundLayer;
		}
	}

	size_t getNeighborBaseEdge(size_t firstNode,size_t neighborIndex)const{
		if(isOutFlow){
			return baseGraph.edgeFromVertex(firstNode,neighborIndex);
		}
		else{
			return baseGraph.edgeToVertex(firstNode,neighborIndex);
		}
	}

	size_t getNeighborBaseVertexRev(size_t firstNode,size_t neighborIndex)const{
		if(!isOutFlow){
			return baseGraph.vertexFromVertex(firstNode,neighborIndex);
		}
		else{
			return baseGraph.vertexToVertex(firstNode,neighborIndex);
		}
	}
	size_t getNeighborLiftedEdge(size_t firstNode,size_t neighborIndex){
		if(isOutFlow){
			return liftedGraph.edgeFromVertex(firstNode,neighborIndex);
		}
		else{
			return liftedGraph.edgeToVertex(firstNode,neighborIndex);
		}
	}
	size_t getNeighborLiftedVertex(size_t firstNode,size_t neighborIndex){
		if(isOutFlow){
			return liftedGraph.vertexFromVertex(firstNode,neighborIndex);
		}
		else{
			return liftedGraph.vertexToVertex(firstNode,neighborIndex);
		}
	}

	size_t numberOfNeighborsLifted(size_t nodeIndex){
		if(isOutFlow){
			return liftedGraph.numberOfEdgesFromVertex(nodeIndex);
		}
		else{
			return liftedGraph.numberOfEdgesToVertex(nodeIndex);
		}
	}

	bool reachable(size_t firstVertex,size_t secondVertex){
		if(isOutFlow){
			return ldpInstance.isReachable(firstVertex,secondVertex);
		}
		else{
			return ldpInstance.isReachable(secondVertex,firstVertex);
		}
	}


	std::pair<bool,size_t> findEdgeBase(size_t firstNode,size_t secondNode){
		if(isOutFlow){
			return baseGraph.findEdge(firstNode,secondNode);
		}
		else{
			return baseGraph.findEdge(secondNode,firstNode);
		}

	}

	size_t getVertexToReach()const{
		if(isOutFlow){
			return ldpInstance.getTerminalNode();
		}
		else{
			return ldpInstance.getSourceNode();
		}
	}


	mutable std::size_t optimalSolutionBase;
	mutable std::unordered_set<size_t> optimalSolutionLifted;

	std::size_t minLayer;
	std::size_t maxLayer;

	const bool isOutFlow;

	const andres::graph::Digraph<>& baseGraph;
	const andres::graph::Digraph<>& liftedGraph;
	const LDP_INSTANCE& ldpInstance;


	std::unordered_map<size_t,double> baseCosts;
	std::unordered_map<size_t,double> liftedCosts;
	mutable std::unordered_map<size_t,double> solutionCosts;
//	mutable std::unordered_map<size_t,size_t> indexStructure;
	//mutable std::unordered_map<size_t,std::unordered_set<size_t>> indexStructure;

	//mutable std::unordered_map<size_t,double> valuesStructure;  //For DFS procedure

	mutable bool optLiftedUpToDate;
	mutable bool optBaseUpToDate;

	mutable double optValue;

	mutable StrForUpdateValues strForUpdateValues;



	//	 std::pair<bool,size_t> findEdgeLifted(size_t firstNode,size_t secondNode){
	//		 if(isOutFlow){
	//			 return liftedGraph.findEdge(firstNode,secondNode);
	//		 }
	//		 else{
	//			 return liftedGraph.findEdge(secondNode,firstNode);
	//		 }
	//	 }



};

template<class LDP_INSTANCE>
inline  ldp_single_node_cut_factor<LDP_INSTANCE>::ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut):
baseGraph(ldpInst.getGraph()),
liftedGraph(ldpInst.getGraphLifted()),
nodeID(nID),
ldpInstance(ldpInst),
isOutFlow(isOut),
nodeNotActive(nID),strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID)
{
	primalBase_=nodeNotActive;  //corresponds to no edge active
	optimalSolutionBase=nodeNotActive;

	if(isOutFlow){
		minLayer=ldpInst.getGroupIndex(nodeID);
		maxLayer=minLayer+ldpInst.getGapLifted(); //some method that returns max time gap lifted
	}
	else{
		maxLayer=ldpInst.getGroupIndex(nodeID);
		minLayer=std::max(0,int(maxLayer)-int(ldpInst.getGapLifted()));
	}

	//baseCosts=std::unordered_map<size_t,double>();
	initBaseCosts(0);
	//liftedCosts=std::unordered_map<size_t,double>();
	initLiftedCosts(0);
	//solutionCosts=std::unordered_map<size_t,double>();
	solutionCosts[nodeNotActive]=0;
	baseCosts[nodeNotActive]=0;
	optLiftedUpToDate=false;
	optBaseUpToDate=false;

	optValue=0;

}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateOptimal() const{
	if(!optLiftedUpToDate) updateValues();
	else if(!optBaseUpToDate){
		optimalSolutionBase=nodeNotActive;
		double minValue=solutionCosts[nodeNotActive];
		for (auto it=solutionCosts.begin();it!=solutionCosts.end();it++) {
			double value=it->second;
			if(value<minValue){
				minValue=value;
				optimalSolutionBase=it->first;
			}
		}
		optBaseUpToDate=true;
	}
	optValue=solutionCosts.at(optimalSolutionBase);
}



template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setPrimalLifted(std::unordered_set<size_t>& verticesOfActiveEdges) {
	primalLifted_=verticesOfActiveEdges;
}




template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::EvaluatePrimal() const{
	double value=0;
	value+=baseCosts.at(primalBase_);
	for(size_t node:primalLifted_){
		value+=liftedCosts.at(node);
	}
	return value;
}




template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setBaseEdgeActive(size_t vertex){
	assert(vertex!=nodeID&&baseCosts.count(vertex)>0);
	primalBase_=vertex;

}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setNoBaseEdgeActive(size_t vertex){

	primalBase_=nodeNotActive;
}



template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getOneBaseEdgeMinMarginal(size_t vertex){
	assert(vertex!=nodeID&&baseCosts.count(vertex)>0);
	updateOptimal();
	if(optimalSolutionBase!=vertex){
		return solutionCosts[vertex]-solutionCosts[optimalSolutionBase];
	}
	else{
		double secondBest=std::numeric_limits<double>::infinity;
		for (auto it=solutionCosts.begin();it!=solutionCosts.end();it++) {
			if(it->first==vertex) continue;
			double value=it->second;
			if(value<secondBest){
				secondBest=value;
			}
		}
		return solutionCosts[vertex]-secondBest;
	}
}


template<class LDP_INSTANCE>
inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllBaseMinMarginals(){
	updateOptimal();
	std::unordered_map<size_t,double> minMarginals;
	if(optimalSolutionBase==nodeNotActive){
		double value=solutionCosts.at(nodeNotActive);
		for(auto pair:solutionCosts){
			if(pair.first!=nodeNotActive){
				minMarginals[pair.first]=pair.second-value;
			}
		}
	}
	else{
		double secondBest=std::numeric_limits<double>::infinity();
		double optValue=solutionCosts.at(optimalSolutionBase);
		for (auto it=solutionCosts.begin();it!=solutionCosts.end();it++) {
			if(it->first==optimalSolutionBase) continue;
			double value=it->second;
			if(value<secondBest){
				secondBest=value;
			}
		}

		for (auto it=solutionCosts.begin();it!=solutionCosts.end();it++) {
			double value=it->second;
			if(it->first!=nodeNotActive){
				minMarginals[it->first]=value-secondBest;
			}
		}
	}
	return minMarginals;
}



template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateCostSimple(const double value,const size_t vertexIndex,bool isLifted){//Only cost change
	if(!isLifted){ //update in base edge
		assert(baseCosts.count(vertexIndex)>0);
		baseCosts[vertexIndex]+=value;
		solutionCosts[vertexIndex]+=value;
		optBaseUpToDate=false;
	}
	else{ //update in lifted edge
		assert(liftedCosts.count(vertexIndex)>0);
		liftedCosts[vertexIndex]+=value;
		//valuesStructure[vertexIndex]+=value;
		optLiftedUpToDate=false;
		optBaseUpToDate=false;
	}
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::addLiftedEdge(size_t node,double cost){
	assert(reachable(nodeID,node));
	liftedCosts[node]=cost;
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::initLiftedCosts(double fractionLifted){
	if(fractionLifted==0){
		for (int i = 0; i < numberOfNeighborsLifted(nodeID); ++i) {
			size_t neighborID=getNeighborLiftedVertex(nodeID,i);
			liftedCosts[neighborID]=0;
		}
	}
	else{
		for (int i = 0; i < numberOfNeighborsLifted(nodeID); ++i) {
			size_t edgeID=getNeighborLiftedEdge(nodeID,i);
			size_t neighborID=getNeighborLiftedVertex(nodeID,i);
			double cost=ldpInstance.getLiftedEdgeScore(edgeID);
			liftedCosts[neighborID]=fractionLifted*cost;
		}
	}
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::initBaseCosts(double fractionBase){
	if(fractionBase==0){
		for (int i = 0; i < numberOfNeighborsBase(nodeID); ++i) {
			size_t neighborID=getNeighborBaseVertex(nodeID,i);
			baseCosts[neighborID]=0;
			solutionCosts[neighborID]=0;
		}
	}
	else{
		for (int i = 0; i < numberOfNeighborsBase(nodeID); ++i) {
			size_t edgeID=getNeighborBaseEdge(nodeID,i);
			size_t neighborID=getNeighborBaseVertex(nodeID,i);
			double cost=ldpInstance.getEdgeScore(edgeID);
			baseCosts[neighborID]=fractionBase*cost;
		}
	}
	//	for (int i = 0; i < numberOfNeighborsLifted(nodeID); ++i) {
	//		size_t edgeID=getNeighborLiftedEdge(nodeID,i);
	//		size_t neighborID=getNeighborLiftedVertex(nodeID,i);
	//		double cost=ldpInstance.getLiftedEdgeScore(edgeID);
	//		liftedCosts[neighborID]=fractionLifted*cost;
	//	}

}


template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::LowerBound() const{//TODO store info about how valuesStructures changed. At least max time layer of changed lifted edge
	updateOptimal();
	return solutionCosts.at(optimalSolutionBase);
}



template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateValues() const{

	//StrForUpdateValues strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID);
	strForUpdateValues.indexStructure.clear();
	strForUpdateValues.solutionCosts.clear();
	strForUpdateValues.valuesStructure.clear();

	updateValues(strForUpdateValues);

	//std::cout<<"values updated"<<std::endl;
	optimalSolutionBase=strForUpdateValues.indexStructure[nodeID];

	size_t vertexInOptimalPath=strForUpdateValues.indexStructure[nodeID];
	optimalSolutionLifted.clear();
	bool hasOptDescendant=vertexInOptimalPath!=nodeNotActive;
	while(hasOptDescendant){
		//std::cout<<"opt vertex "<<vertexInOptimalPath<<std::endl;
		if(liftedCosts.count(vertexInOptimalPath)>0){
			optimalSolutionLifted.insert(vertexInOptimalPath);
		}
		auto it=strForUpdateValues.indexStructure.find(vertexInOptimalPath);
		hasOptDescendant=it!=strForUpdateValues.indexStructure.end();
		if(hasOptDescendant){
			vertexInOptimalPath=it->second;
		}
	}


	optLiftedUpToDate=true;
	optBaseUpToDate=true;
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateValues(StrForUpdateValues& myStr,size_t vertexToIgnore) const{

	//std::cout<<"update values in node "<<nodeID<<std::endl;
	std::unordered_set<size_t> closedVertices;


	bool lastLayerSet=false;
	size_t lastLayer=0;
	if(liftedCosts.count(vertexToIgnore)>0){
		lastLayerSet=true;
		lastLayer=ldpInstance.getGroupIndex(vertexToIgnore);
	}


	std::stack<size_t> nodeStack;
	nodeStack.push(nodeID);

	while(!nodeStack.empty()){
		size_t currentNode=nodeStack.top();
		//std::cout<<"current node "<<currentNode<<std::endl;
		bool descClosed=true;
		double minValue=0;
		//std::unordered_set<size_t> minValueIndices;
		size_t minValueIndex=getVertexToReach();

		for (int i = 0; i < numberOfNeighborsBase(currentNode); ++i) {
			size_t desc=getNeighborBaseVertex(currentNode,i);
			if(desc==vertexToIgnore||desc==getVertexToReach()) continue;
			if(isInThisFactorRange(desc)&&myStr.useVertex(desc)){
				if(closedVertices.count(desc)>0||(lastLayerSet&&!isInGivenRange(desc,lastLayer))){  //descendant closed
					if(descClosed){
						auto it=myStr.valuesStructure.find(desc);
						if(it!=myStr.valuesStructure.end()){
							if(minValue>it->second){
								minValue=it->second;
								minValueIndex=desc;
							}
						}
					}
				}
				else {  //descendant not closed
					nodeStack.push(desc);
					descClosed=false;

				}
			}

		}
		if(descClosed){ //Close node if all descendants are closed
		//	std::cout<<"desc closed"<<std::endl;
			if(currentNode==nodeID){  //all nodes closed, compute solution values
				double bestValue=0;
				size_t bestVertex=nodeNotActive;
				for (auto it=myStr.baseCosts.begin();it!=myStr.baseCosts.end();++it) {
				//	std::cout<<"base edge "<<it->first<<std::endl;
					size_t vertex=it->first;
					double baseCost=it->second;
					double valueToAdd=0;
					auto vsIt=myStr.valuesStructure.find(vertex);
					if(vsIt!=myStr.valuesStructure.end()){
						valueToAdd=vsIt->second;
					}
					else{
						const auto lcIt=myStr.liftedCosts.find(vertex);
						if(lcIt!=myStr.liftedCosts.end()){
							valueToAdd=lcIt->second;
						}
					}
					double value=baseCost+valueToAdd;
					myStr.solutionCosts[vertex]=value;
					if(value<bestValue){
						bestValue=value;
						bestVertex=vertex;
					}
					//std::cout<<"end for"<<std::endl;
				}
				//std::cout<<"after for"<<std::endl;
				//optimalSolution=bestVertex;
				myStr.indexStructure[nodeID]=bestVertex;
				myStr.optValue=bestValue;
			}
			else{

				double valueToStore=minValue;
				const auto liftedIt=myStr.liftedCosts.find(currentNode);
				if(liftedIt!=myStr.liftedCosts.end()){
					valueToStore+=liftedIt->second;  //add lifted edge cost if the edge exists
				}
				const auto baseIt=baseCosts.find(currentNode);
				if(valueToStore<=0||(baseIt!=myStr.baseCosts.end()&&minValue<0)){  //store only negative values or values needed to correct solutionCosts
					myStr.valuesStructure[currentNode]=valueToStore;
					myStr.indexStructure[currentNode]=minValueIndex;
				}
//				else{
//					myStr.valuesStructure.erase(currentNode);
//					myStr.indexStructure.erase(currentNode);
//				}
				closedVertices.insert(currentNode); //marking the node as closed.
			}
			nodeStack.pop();
			//std::cout<<"pop from node stack"<<std::endl;
		}
	}



}




class ldp_mcf_single_node_cut_message
{
	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamLeft(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
	{
		r.updateCostSimple(msg,msg_dim);
		//Only base edges are updated if message comes from mcf, the dfs procedure not needed

	}

	template<typename MCF_FACTOR>
	void RepamRight(MCF_FACTOR& r, const double msg, const std::size_t msg_dim) const
	{
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
	{
	}

	template<typename MCF_FACTOR, typename MSG_ARRAY>
	static void SendMessagesToRight(const MCF_FACTOR& leftRepam, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
	{
	}
};


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::findAllOptimal(std::unordered_set<size_t>& isNotZeroInOpt,std::unordered_set<size_t>& isOneInOpt,std::unordered_map<size_t,std::unordered_set<size_t>>& candidateGraph,const StrForUpdateValues& strForUpdateValues){
	std::stack<size_t> myStack;
	std::list<size_t> parentStack;

	myStack.push(nodeID);
	parentStack.push_back(nodeID);

	//std::unordered_set<size_t> descendants;
	double bestValue=strForUpdateValues.optValue;
	for(auto pair:strForUpdateValues.solutionCosts){
		size_t desc=pair.first;
		double value=pair.second;
		if(value==bestValue){
			//descendants.insert(desc);
			if(desc!=nodeNotActive&&desc!=getVertexToReach()){
				myStack.push(desc);
			}
			else{
				isNotZeroInOpt.clear();//TODO find out if this is necessary
			}
		}
	}
//	for(size_t desc:strForUpdateValues.solutionCosts){
//		myStack.push(desc);
//	}
	while(!myStack.empty()){
		isOneInOpt.insert(myStack.top());
		if(*(parentStack.rbegin())==myStack.top()){
			parentStack.pop_back();
			myStack.pop();
		}
		else{
			size_t currentVertex=myStack.top();

			double bestValue=0;
			const auto descIt=strForUpdateValues.indexStructure.find(currentVertex);
			if(descIt!=strForUpdateValues.indexStructure.end()){
				size_t desc=descIt->second;
				if(desc!=getVertexToReach()){
					bestValue=strForUpdateValues.valuesStructure.at(desc);
					myStack.push(desc);
				}
			}

			std::unordered_set<size_t>& candidates=candidateGraph[currentVertex];
			if(!candidates.empty()){
				parentStack.push_back(currentVertex);
				for(size_t desc:candidates){
					if(desc==descIt->second) continue;

					double value=strForUpdateValues.valuesStructure.at(desc);
					if(value==bestValue){
						myStack.push(desc);
					}

				}
			}
			if(bestValue==0){  //The optimal path can be terminated in the current vertex
				if(!isNotZeroInOpt.empty()){
					std::unordered_set<size_t> keepOpen;
					keepOpen.insert(currentVertex);
					for(size_t optVertex:parentStack){
						if(isNotZeroInOpt.count(optVertex)>0){
							keepOpen.insert(optVertex);
						}
					}
					isNotZeroInOpt=keepOpen;
				}
				myStack.pop();
			}
		}
	}
}



template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::oneLiftedMinMarginal(size_t vertexOfLiftedEdge){
	updateOptimal();

	std::unordered_set<size_t> isOneInOpt(optimalSolutionLifted);
	std::unordered_set<size_t> isNotZeroInOpt(optimalSolutionLifted);
	std::unordered_map<size_t,std::unordered_set<size_t>> candidateGraph=createCandidateGraph(strForUpdateValues);
	findAllOptimal(isNotZeroInOpt,isOneInOpt,candidateGraph,strForUpdateValues);

	if(isOneInOpt.count(vertexOfLiftedEdge)>0){
		if(isNotZeroInOpt.count(vertexOfLiftedEdge)>0){
			std::unordered_map<size_t,double> localSolutionCosts;
			StrForUpdateValues myStr(baseCosts,liftedCosts,localSolutionCosts,nodeID);
			updateValues(myStr,vertexOfLiftedEdge);
			//size_t restrictedOptimalSolution=myStr.indexStructure[nodeID];
			double restrictedOptValue=myStr.optValue;
			return optValue-restrictedOptValue;
		}
		else{
			return 0;
		}
	}
	else{
		std::unordered_map<size_t,double> message=bottomUpUpdate(strForUpdateValues,vertexOfLiftedEdge);
		auto it =message.begin();

		return it->second;
	}

}


template<class LDP_INSTANCE>
inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::bottomUpUpdate(const StrForUpdateValues& myStr,size_t vertex,std::unordered_set<size_t>* pClosedVert,std::unordered_map<size_t,double>* pBUValuesStr){
	bool onlyOne=pClosedVert==0;

	std::unordered_map<size_t,double> messages;
	if(onlyOne){
		pClosedVert=new std::unordered_set<size_t> ();
		pBUValuesStr=new std::unordered_map<size_t,double>();
	}
	std::unordered_set<size_t> closedVertices=*pClosedVert;
	std::unordered_map<size_t,double> buValuesStr=*pBUValuesStr;

	std::stack<size_t> myStack;
	myStack.push(vertex);
	while(!myStack.empty()){
		size_t currentVertex=myStack.top();
		bool predClosed=true;
		for (int i = 0; i < numberOfNeighborsBaseRev(currentVertex); ++i) {
			size_t pred=getNeighborBaseVertexRev(currentVertex,i);
			if(pred==nodeID) continue;
			if(reachable(nodeID,pred)&&closedVertices.count(pred)==0){
				predClosed=false;
				myStack.push(pred);
			}
		}
		if(predClosed){
			double bestValue=std::numeric_limits<double>::infinity();
			auto baseIt=myStr.baseCosts.find(currentVertex);
			if(baseIt!=myStr.baseCosts.end()){
				bestValue=baseIt->second;
			}
			for (int i = 0; i < numberOfNeighborsBaseRev(currentVertex); ++i) {
				size_t pred=getNeighborBaseVertexRev(currentVertex,i);
				if(pred==nodeID) continue;
				auto valuesIt=buValuesStr.find(pred);
				if(valuesIt!=buValuesStr.end()){
					double value=valuesIt->second;
					if(value<bestValue){
						bestValue=value;
					}
				}
			}
			auto liftedIt=myStr.liftedCosts.find(currentVertex);
			if(liftedIt!=myStr.liftedCosts.end()){
				if(onlyOne&&currentVertex!=vertex){
					bestValue+=liftedIt->second;
				}
				else{
					double topDownValue=0;
					auto valuesIt=myStr.valuesStructure.find(currentVertex);
					if(valuesIt!=myStr.valuesStructure.end()){
						topDownValue=valuesIt->second;
					}
					double restrictedOpt=topDownValue+bestValue;
					bestValue+=liftedIt->second;

					double delta=restrictedOpt-myStr.optValue;
					messages[currentVertex]=delta;
					bestValue-=delta;
					//myStr.liftedCosts-=delta;

				}
			}
			closedVertices.insert(currentVertex);
			buValuesStr[currentVertex]=bestValue;
			myStack.pop();

		}

	}
	if(onlyOne){
		delete(pClosedVert);
		delete(pBUValuesStr);
		pClosedVert=0;
		pBUValuesStr=0;
	}
	return messages;


}

template<class LDP_INSTANCE>
inline std::unordered_map<size_t,std::unordered_set<size_t>> ldp_single_node_cut_factor<LDP_INSTANCE>::createCandidateGraph(const StrForUpdateValues& myStr){
	std::unordered_map<size_t,std::unordered_set<size_t>> candidateGraph;
		for(auto it=myStr.valuesStructure.begin();it!=myStr.valuesStructure.end();++it){
			size_t vertex=it->first;
			//strForUpdateValues.relevantVertices.insert(vertex);
			for (int i = 0; i < numberOfNeighborsBase(vertex); ++i) {
				size_t neighbor=getNeighborBaseVertex(vertex,i);
				if(myStr.valuesStructure.count(neighbor)>0){
					candidateGraph[vertex].insert(neighbor);
				}
			}
		}

		return candidateGraph;
}

template<class LDP_INSTANCE>
inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllLiftedMinMarginals(){

	//std::unordered_map<size_t,double> baseMessages=getAllBaseMinMarginals();
	updateOptimal();
	std::unordered_map<size_t,double> liftedMessages;

	std::unordered_set<size_t> isNotZeroInOpt;
	std::unordered_set<size_t> isOneInOpt;

	std::unordered_map<size_t,double> localSolutionCosts=solutionCosts;
	std::unordered_map<size_t,double> localLiftedCosts=liftedCosts;
	std::unordered_map<size_t,double> localBaseCosts=baseCosts;
//	for(auto pair:baseMessages){
//		localBaseCosts[pair.first]-=pair.second;
//		localSolutionCosts[pair.first]-=pair.second;
//	}

	StrForUpdateValues myStr(localBaseCosts,localLiftedCosts,localSolutionCosts,nodeID);

	double currentOptValue=myStr.optValue;

	//updateValues();
	std::list<size_t> openVertices;

	size_t vertexInOptimalPath=myStr.indexStructure[nodeID];
	//openVertices.push_back(nodeID);


	//size_t vertexInOptimalPath=strForUpdateValues.indexStructure[nodeID];
	bool hasOptDescendant=vertexInOptimalPath!=nodeNotActive;
	while(hasOptDescendant){
		isOneInOpt.insert(vertexInOptimalPath);  //Probably not necessary?
		//isNotZeroInOpt.insert(vertexInOptimalPath);
		if(liftedCosts.count(vertexInOptimalPath)>0){
			openVertices.push_back(vertexInOptimalPath);
			isNotZeroInOpt.insert(vertexInOptimalPath);
		}
		auto it=myStr.indexStructure.find(vertexInOptimalPath);
		hasOptDescendant=it!=myStr.indexStructure.end();
		if(hasOptDescendant){
			vertexInOptimalPath=it->second;
		}
	}


	std::unordered_map<size_t,std::unordered_set<size_t>> candidateGraph=createCandidateGraph(myStr);

	findAllOptimal(isNotZeroInOpt,isOneInOpt,candidateGraph,myStr);

	myStr.setUseAllVertices(false);


	for (auto listIt=openVertices.begin();listIt!=openVertices.end();listIt++) {
		size_t vertexToClose=*listIt;
		if(isNotZeroInOpt.count(vertexToClose)>0){

			updateValues(myStr,vertexToClose);
			double newOpt=myStr.optValue;
			findAllOptimal(isNotZeroInOpt,isOneInOpt,candidateGraph,myStr);  //valuesStructure and indexStructure can be used global
			double delta=currentOptValue-newOpt;
			localLiftedCosts[vertexToClose]-=delta;
			liftedMessages[vertexToClose]+=delta;
			currentOptValue=newOpt;
			isNotZeroInOpt.erase(vertexToClose);
		}
	}
	updateValues(myStr);
	std::unordered_map<size_t,double> buValuesStructure;
	std::unordered_set<size_t> closedVertices;
	for(size_t optVertex:isOneInOpt){
		buValuesStructure[optVertex]=currentOptValue-myStr.valuesStructure[optVertex]+localLiftedCosts[optVertex];
		closedVertices.insert(optVertex);
	}
	for(auto pair:liftedCosts){
		if(closedVertices.count(pair.first)==0){
			std::unordered_map<size_t,double> newMessages=bottomUpUpdate(myStr,pair.first,&closedVertices,&buValuesStructure);
			liftedMessages.insert(newMessages.begin(),newMessages.end());
		}
	}


	return liftedMessages;
}








class ldp_mcf_single_node_cut_base_edge_message
{
	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamLeft(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
	{
		assert(msg_dim == 0);
		r.updateCostSimple(msg,right_node);
	}

	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamRight(SINGLE_NODE_CUT_FACTOR& l, const double msg, const std::size_t msg_dim) const
	{
		assert(msg_dim == 0);
		l.updateCostSimple(msg,left_node);
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
	{
	}

	template<typename MCF_FACTOR, typename MSG_ARRAY>
	static void SendMessagesToRight(const MCF_FACTOR& leftRepam, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
	{
	}

	template<typename SINGLE_NODE_CUT_FACTOR>
	bool check_primal_consistency(const SINGLE_NODE_CUT_FACTOR& l, const SINGLE_NODE_CUT_FACTOR& r) const
	{

	}

private:
	std::size_t left_node;
	std::size_t right_node;
};






}
