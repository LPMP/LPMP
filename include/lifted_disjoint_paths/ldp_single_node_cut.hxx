#pragma once
#include <andres/graph/digraph.hxx>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <array>


namespace LPMP {

struct baseAndLiftedMessages{
	std::unordered_map<size_t,double> baseMessages;
	std::unordered_map<size_t,double> liftedMessages;
};

struct StrForUpdateValues{
	std::unordered_map<size_t,double> solutionCosts;
	std::unordered_map<size_t,double> valuesStructure;
	std::unordered_map<size_t,size_t> indexStructure;
	const std::unordered_map<size_t,double>& baseCosts;
	const std::unordered_map<size_t,double>& liftedCosts;
	std::unordered_set<size_t> relevantVertices;
	bool useAllVertices;
	//size_t optimalSolution;
	StrForUpdateValues(const std::unordered_map<size_t,double>& bCosts,const std::unordered_map<size_t,double>& lCosts):
	baseCosts(bCosts),
	liftedCosts(lCosts)
	{
		useAllVertices=true;
		solutionCosts=baseCosts;
		//optimalSolution=nodeID;
	}
	bool useVertex(size_t vertex){
		if(useAllVertices) return true;
		return relevantVertices.count(vertex)>0;
	}

	void setUseAllVertices(bool value){
		useAllVertices=value;
	}
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


	template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primalBaseComplete_); }
	template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

	//auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }//TODO change this. This will not work with so many variables
	double tmp_to_delete_val;
	auto export_variables() { return std::tie(tmp_to_delete_val); } //?What comes here?

	void updateCostSimple(const double value,const size_t vertexIndex,bool isLifted);
	double getBaseEdgeMinMarginal(size_t vertex);
	std::unordered_map<size_t,double> getAllBaseMinMarginals();

	//std::unordered_map<size_t,double> adjustCostsAndSendMessages();



	const std::size_t nodeID;
	const size_t nodeNotActive;
	size_t primalBase_;
	std::unordered_set<size_t> primalLifted_;

	//	void init_primal() {
	//		for(auto pair:baseCosts){
	//			size_t vertex=pair.first;
	//			if(vertex!=nodeID){
	//				auto fe=findEdgeBase(nodeID,vertex);
	//				size_t edgeIndex=fe.second;
	//				primal_[edgeIndex]=0;
	//			}
	//		}
	//	}

//	static void init_primal_vector(andres::graph::Digraph<> graph) {
//		primalBaseComplete_=std::vector<bool>(graph.numberOfEdges(),0);
//	}
//
//	static void init_primal_vector(size_t numberOfEdges) {
//		primalBaseComplete_=std::vector<bool>(numberOfEdges,0);
//	}


	//	double getLiftedEdgeMinMarginal();


	baseAndLiftedMessages adjustCostsAndSendMessagesLifted();

	//std::unordered_map<size_t,std::unordered_set<size_t>> initPredecessorsInIndexStr();



private:
	void updateValues() const;
	void updateValues(StrForUpdateValues& myStr,size_t lastLayer=0) const;
	void updateOptimal() const;
	void findAllOptimal(std::unordered_set<size_t>& isNotZeroInOpt,std::unordered_set<size_t>& isOneInOpt,std::unordered_map<size_t,std::unordered_set<size_t>>& candidateGraph,StrForUpdateValues& strForUpdateValues);

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


	mutable std::unordered_map<size_t,double> baseCosts;
	std::unordered_map<size_t,double> liftedCosts;
	mutable std::unordered_map<size_t,double> solutionCosts;
	mutable std::unordered_map<size_t,size_t> indexStructure;
	//mutable std::unordered_map<size_t,std::unordered_set<size_t>> indexStructure;

	//mutable std::unordered_map<size_t,double> valuesStructure;  //For DFS procedure

	mutable bool vsUpToDate;
	mutable bool optUpToDate;


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
nodeNotActive(nID)
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
	vsUpToDate=false;
	optUpToDate=false;

}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateOptimal() const{
	if(!vsUpToDate) updateValues();
	else if(!optUpToDate){
		optimalSolutionBase=nodeNotActive;
		double minValue=solutionCosts[nodeNotActive];
		for (auto it=solutionCosts.begin();it!=solutionCosts.end();it++) {
			double value=it->second;
			if(value<minValue){
				minValue=value;
				optimalSolutionBase=it->first;
			}
		}
		optUpToDate=true;
	}
}



template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setPrimalLifted(std::unordered_set<size_t>& verticesOfActiveEdges) {
	primalLifted_=verticesOfActiveEdges;
}


//template<class LDP_INSTANCE>
//inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setPrimalLifted() {
//	double value=0;
//	//bool search=true;
//	primalLifted_.clear();
//	if(primalBase_!=nodeNotActive){
//		size_t currentVertex=primalBase_;
//		if(liftedCosts.count(primalBase_)>0){
//			primalLifted_.insert(primalBase_);
//		}
//		bool search=isInThisFactorRange(currentVertex)&&currentVertex!=getVertexToReach();
//		while(search){
//			size_t numberOfNeighbors=numberOfNeighborsBase(currentVertex);
//			for (size_t i = 0; i < numberOfNeighbors; ++i) {
//				size_t edgeId=getNeighborBaseEdge(currentVertex,i);
//				if(primalBaseComplete_[edgeId]==1){
//					currentVertex=getNeighborBaseVertex(currentVertex,i);
//					if(liftedCosts.count(currentVertex)>0){
//						primalLifted_.insert(currentVertex);
//					}
//					break;
//				}
//			}
//			if(!isInThisFactorRange(currentVertex)||currentVertex==getVertexToReach()){
//				search=false;
//			}
//		}
//	}
//}

template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::EvaluatePrimal() const{
	double value=0;
	value+=baseCosts[primalBase_];
	for(size_t node:primalLifted_){
		value+=liftedCosts.at(node);
	}
	return value;
}


//template<class LDP_INSTANCE>
//inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setBaseEdgeLabel(size_t vertex,bool value){
//	if(value){
//		setBaseEdgeActive(vertex);
//	}
//	else{
//		setBaseEdgeInactive(vertex);
//	}
//}

template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setBaseEdgeActive(size_t vertex){
	assert(vertex!=nodeID&&baseCosts.count(vertex)>0);
	primalBase_=vertex;
//	for (int i = 0; i < numberOfNeighborsBase(nodeID); ++i) {
//		size_t edgeIndex=getNeighborBaseEdge(nodeID,i);
//		primalBaseComplete_[edgeIndex]=0;
//	}
//	auto fe=findEdgeBase(nodeID,vertex);
//	size_t edgeIndex=fe.second;
//	primalBaseComplete_[edgeIndex]=1;
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setNoBaseEdgeActive(size_t vertex){
//	for (int i = 0; i < baseCosts.size()-1; ++i) {
//		size_t edgeID=getNeighborBaseEdge(nodeID,i);
//		primalBaseComplete_[edgeID]=0;
//	}
	primalBase_=nodeNotActive;
}



template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getBaseEdgeMinMarginal(size_t vertex){
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
		double optValue=solutionCosts[optimalSolutionBase];
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
//
//	std::array<double,2> twoBestValues;
//	twoBestValues={solutionCosts.at(nodeNotActive),std::numeric_limits<double>::infinity()};
//	optimalSolution=nodeNotActive;
//
//	for (auto it=solutionCosts.begin();it!=solutionCosts.end();it++) {
//		double value=it->second;
//		if(value<twoBestValues[0]){
//			twoBestValues[1]=twoBestValues[0];
//			twoBestValues[0]=it->second;
//			optimalSolution=it->first;
//		}
//		else if(value<twoBestValues[1]){
//			twoBestValues[1]=value;
//		}
//	}
//	for (auto it=solutionCosts.begin();it!=solutionCosts.end();it++) {
//		if(it->first==optimalSolution){
//			minMarginals[optimalSolution]=twoBestValues[0]-twoBestValues[1];
//		}
//		else{
//			minMarginals[it->first]=it->second-twoBestValues[0];
//		}
//	}
//	optUpToDate=true;
//	return minMarginals;
//}



//template<class LDP_INSTANCE>
//inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::adjustCostsAndSendMessages(){
//	updateOptimal();
//	double minValue=solutionCosts[optimalSolution];
//	std::unordered_map<size_t,double> messages;
//	for (auto it=solutionCosts.begin();it!=solutionCosts.end();it++) {
//		if(it->first==optimalSolution) continue;
//		double delta=it->second-minValue;
//		messages[it->first]=delta;
//		baseCosts[it->first]-=delta;
//		it->second=minValue;
//	}
//	return messages;
//}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateCostSimple(const double value,const size_t vertexIndex,bool isLifted){//Only cost change
	if(!isLifted){ //update in base edge
		assert(baseCosts.count(vertexIndex)>0);
		baseCosts[vertexIndex]+=value;
		solutionCosts[vertexIndex]+=value;
		optUpToDate=false;
	}
	else{ //update in lifted edge
		assert(liftedCosts.count(vertexIndex)>0);
		liftedCosts[vertexIndex]+=value;
		//valuesStructure[vertexIndex]+=value;
		vsUpToDate=false;
		optUpToDate=false;
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

	StrForUpdateValues strForUpdateValues(baseCosts,liftedCosts);
	updateValues(strForUpdateValues);

	solutionCosts=strForUpdateValues.solutionCosts;
	optimalSolutionBase=strForUpdateValues.indexStructure[nodeID];

	size_t vertexInOptimalPath=strForUpdateValues.indexStructure[nodeID];
	optimalSolutionLifted.clear();
	bool hasOptDescendant=vertexInOptimalPath!=nodeNotActive;
	while(hasOptDescendant){
		if(liftedCosts.count(vertexInOptimalPath)>0){
			optimalSolutionLifted.insert(vertexInOptimalPath);
		}
		hasOptDescendant=strForUpdateValues.indexStructure.count(vertexInOptimalPath)>0;
		if(hasOptDescendant){
			size_t newVertex=strForUpdateValues.indexStructure[vertexInOptimalPath];
			vertexInOptimalPath=newVertex;
		}
	}


	vsUpToDate=true;
	optUpToDate=true;
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateValues(StrForUpdateValues& myStr,size_t lastLayer) const{

	std::unordered_set<size_t> closedVertices;

	bool lastLayerSet=lastLayer!=0;

	std::stack<size_t> nodeStack;
	nodeStack.push(nodeID);

	while(!nodeStack.empty()){
		size_t currentNode=nodeStack.top();
		bool descClosed=true;
		double minValue=0;
		std::unordered_set<size_t> minValueIndices;
		size_t minValueIndex=nodeNotActive;

		for (int i = 0; i < numberOfNeighborsBase(currentNode); ++i) {
			size_t desc=getNeighborBaseVertex(currentNode,i);
			if((!lastLayerSet&&isInThisFactorRange(desc))||(lastLayerSet&&isInGivenRange(desc,lastLayer))){
				if(myStr.useVertex(desc)){
					if(closedVertices.count(desc)>0){  //descendant closed
						if(descClosed&&myStr.valuesStructure.count(desc)>0){
							if(minValue>=myStr.valuesStructure[desc]&&myStr.valuesStructure[desc]<0){
								minValue=myStr.valuesStructure[desc];
								minValueIndex=desc;
							}
						}
					}
					else{  //descendant not closed
						nodeStack.push(desc);
						descClosed=false;

					}
				}
			}

		}
		if(descClosed){ //Close node if all descendants are closed
			if(currentNode==nodeID){  //all nodes closed, compute solution values
				double bestValue=0;
				size_t bestVertex=nodeNotActive;
				for (auto it=myStr.baseCosts.begin();it!=myStr.baseCosts.end();it++) {
					size_t vertex=it->first;
					double baseCost=it->second;
					double valueToAdd=0;
					if(myStr.valuesStructure.count(vertex)>0){
						valueToAdd=myStr.valuesStructure[vertex];
					}
					else{
						if(myStr.liftedCosts.count(vertex)>0){
							valueToAdd=myStr.liftedCosts.at(vertex);
						}
					}
					myStr.solutionCosts[vertex]=baseCost+valueToAdd;
					if(myStr.solutionCosts[vertex]<bestValue){
						bestValue=myStr.solutionCosts.at(vertex);
						bestVertex=vertex;
					}
				}
				//optimalSolution=bestVertex;
				myStr.indexStructure[nodeID]=bestVertex;
			}
			else{

				double valueToStore=minValue;
				if(myStr.liftedCosts.count(currentNode)>0){
					valueToStore+=myStr.liftedCosts.at(currentNode);  //add lifted edge cost if the edge exists
				}
				if(valueToStore<0||(myStr.baseCosts.count(currentNode)>0&&minValue<0)){  //store only negative values or values needed to correct solutionCosts
					myStr.valuesStructure[currentNode]=valueToStore;
					myStr.indexStructure[currentNode]=minValueIndex;
				}
				else{
					myStr.valuesStructure.erase(currentNode);
				}
				closedVertices.insert(currentNode); //marking the node as closed.
			}
			nodeStack.pop();
		}
	}



}



//template<class LDP_INSTANCE>
//inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateCostFull(const double value,const size_t index){//Includes DFS update
//	size_t myIndex=decodeIndex(index);
//	if(myIndex<numberOfEdges){ //update in base edge
//		baseCosts[myIndex]+=value;
//		solutionCosts[myIndex]+=value;
//		if(myIndex==optimalSolution){
//			if(value>0){
//				for (int i = 0; i < numberOfEdges; ++i) {
//					if(solutionCosts[i]<solutionCosts[optimalSolution]){
//						optimalSolution=i;
//						primal_=i;
//					}
//				}
//			}
//		}
//		else if(solutionCosts[myIndex]<solutionCosts[optimalSolution]){
//			optimalSolution=myIndex;
//			primal_=myIndex; //maybe not
//
//		}
//
//	}
//	else{ //update in lifted edge
//		liftedCosts[myIndex-numberOfEdges]+=value;
//		updateValues(myIndex-numberOfEdges);
//	}
//}





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
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::findAllOptimal(std::unordered_set<size_t>& isNotZeroInOpt,std::unordered_set<size_t>& isOneInOpt,std::unordered_map<size_t,std::unordered_set<size_t>>& candidateGraph,StrForUpdateValues& strForUpdateValues){
	std::stack<size_t> myStack;
	std::list<size_t> parentStack;

	myStack.push(nodeID);
	parentStack.push_back(nodeID);

	std::unordered_set<size_t> descendants;
	double bestValue=0;
	for(auto pair:strForUpdateValues.solutionCosts){
		size_t desc=pair.first;
		double value=strForUpdateValues.valuesStructure[desc];
		if(value<bestValue){
			descendants.clear();
			descendants.insert(desc);
			bestValue=value;
		}
		else if(value==bestValue){
			descendants.insert(desc);
		}
	}
	for(size_t desc:strForUpdateValues.solutionCosts){
		myStack.push(desc);
	}
	while(!myStack.empty()){
		isOneInOpt.insert(myStack.top());
		if(*(parentStack.rbegin())==myStack.top()){
			parentStack.pop_back();
			myStack.pop();
		}
		else{
			size_t currentVertex=myStack.top();
			double bestValue=0;
			std::unordered_set<size_t> descendants;
			parentStack.push_back(currentVertex);
			std::unordered_set<size_t>& candidates=candidateGraph[currentVertex];
			if(!candidates.empty()){
				for(size_t desc:candidates){
					double value=strForUpdateValues.valuesStructure[desc];
					if(value<bestValue){
						descendants.clear();
						descendants.insert(desc);
						bestValue=value;
					}
					else if(value==bestValue){
						descendants.insert(desc);
					}
				}
				for(size_t desc:descendants){
					myStack.push(desc);
				}
			}
			else{
				std::unordered_set<size_t> keepOpen;
				keepOpen.insert(currentVertex);
				for(size_t optVertex:parentStack){
					if(isNotZeroInOpt.count(optVertex)>0){
						keepOpen.insert(optVertex);
					}
				}
				isNotZeroInOpt=keepOpen;
			}
		}
	}
}



template<class LDP_INSTANCE>
inline baseAndLiftedMessages ldp_single_node_cut_factor<LDP_INSTANCE>::adjustCostsAndSendMessagesLifted(){

	std::unordered_map<size_t,double> baseMessages=getAllBaseMinMarginals();
	std::unordered_map<size_t,double> liftedMessages;

	std::unordered_set<size_t> isNotZeroInOpt;
	std::unordered_set<size_t> isOneInOpt;

	std::unordered_map<size_t,double> localSolutionCosts=solutionCosts;
	std::unordered_map<size_t,double> localLiftedCosts=liftedCosts;
	std::unordered_map<size_t,double> localBaseCosts=baseCosts;
	for(auto pair:baseMessages){
		localBaseCosts[pair.first]-=pair.second;
		localSolutionCosts[pair.first]-=pair.second;
	}

	StrForUpdateValues strForUpdateValues(localBaseCosts,localLiftedCosts);
	strForUpdateValues.solutionCosts=localSolutionCosts;

	double currentOptValue=strForUpdateValues.solutionCosts[optimalSolutionBase];

	//updateValues();
	std::list<size_t> openVertices;

	size_t vertexInOptimalPath=strForUpdateValues.indexStructure[nodeID];
	//openVertices.push_back(nodeID);


	size_t vertexInOptimalPath=strForUpdateValues.indexStructure[nodeID];
	bool hasOptDescendant=vertexInOptimalPath!=nodeNotActive;
	while(hasOptDescendant){
		isOneInOpt.insert(vertexInOptimalPath);
		isNotZeroInOpt.insert(vertexInOptimalPath);
		if(liftedCosts.count(vertexInOptimalPath)>0){
			openVertices.push_back(vertexInOptimalPath);
		}
		hasOptDescendant=strForUpdateValues.indexStructure.count(vertexInOptimalPath)>0;
		if(hasOptDescendant){
			size_t newVertex=strForUpdateValues.indexStructure[vertexInOptimalPath];
			vertexInOptimalPath=newVertex;
		}
	}


	std::unordered_map<size_t,std::unordered_set<size_t>> candidateGraph;
	for(auto it=strForUpdateValues.valuesStructure.begin();it!=strForUpdateValues.valuesStructure.end();it++){
		size_t vertex=it->first;
		strForUpdateValues.relevantVertices.insert(vertex);
		for (int i = 0; i < numberOfNeighborsBase(vertex); ++i) {
			size_t neighbor=getNeighborBaseVertex(vertex,i);
			if(strForUpdateValues.valuesStructure.count(neighbor)>0){
				candidateGraph[vertex].insert(neighbor);
			}
		}
	}

	findAllOptimal(isNotZeroInOpt,isOneInOpt,candidateGraph,strForUpdateValues);

	strForUpdateValues.setUseAllVertices(false);


	for (int i = 0; i < openVertices.size(); ++i) {
		size_t vertexToClose=openVertices[i];
		if(isNotZeroInOpt.count(vertexToClose)>0){
			strForUpdateValues.relevantVertices.erase(vertexToClose);
			updateValues(strForUpdateValues,ldpInstance.getGroupIndex(vertexToClose));
			double newOpt=strForUpdateValues.solutionCosts[nodeID];
			findAllOptimal(isNotZeroInOpt,isOneInOpt,candidateGraph,strForUpdateValues);  //valuesStructure and indexStructure can be used global
			double delta=currentOptValue-newOpt;
			localLiftedCosts[vertexToClose]-=delta;
			liftedMessages[vertexToClose]+=delta;
			currentOptValue-=delta;
			strForUpdateValues.relevantVertices.insert(vertexToClose);
			isNotZeroInOpt.erase(vertexToClose);
		}
	}




	baseAndLiftedMessages messages={baseMessages,liftedMessages};
	return messages;
}













//template<class LDP_INSTANCE>
//inline baseAndLiftedMessages ldp_single_node_cut_factor<LDP_INSTANCE>::adjustCostsAndSendMessagesLifted(){
//	std::unordered_map<size_t,double> baseMessages=adjustCostsAndSendMessages();
//	std::unordered_map<size_t,double> liftedMessages;
//
//	std::unordered_set<size_t> isNotZeroInOpt;
//	std::unordered_set<size_t> isOneInOpt;
//
//	updateValues();
//	std::list<size_t> openVertices;
//
//	size_t vertexInOptimalPath=nodeID;
//	openVertices.push_back(nodeID);
//	std::unordered_map<size_t,size_t> openBackward;
//	std::unordered_map<size_t,size_t> openForward;
//
//
//	while(indexStructure.count(vertexInOptimalPath)>0){
//		isOneInOpt.insert(vertexInOptimalPath);
//		isNotZeroInOpt.insert(vertexInOptimalPath);
//		openBackward[vertexInOptimalPath]=*openVertices.rbegin();
//		if(liftedCosts.count(vertexInOptimalPath)>0){
//			openVertices.push_back(vertexInOptimalPath);
//		}
//		size_t newVertex=indexStructure[vertexInOptimalPath];
//		vertexInOptimalPath=newVertex;
//	}
//
//	std::unordered_map<size_t,std::unordered_set<size_t>> candidateGraph;
//	for(auto it=valuesStructure.begin();it!=valuesStructure.end();it++){
//		size_t vertex=it->first;
//		for (int i = 0; i < numberOfNeighborsBase(vertex); ++i) {
//			size_t neighbor=getNeighborBaseVertex(vertex,i);
//			if(indexStructure.count(neighbor)>0){
//				candidateGraph[vertex].insert(neighbor);
//			}
//		}
//	}
//
//	size_t boundLayer=maxLayer;
//	if(!isOutFlow){
//		boundLayer=minLayer;
//	}
//	std::unordered_set<size_t> isProcessed;
//	while(openVertices.size()>0){
//		std::stack<size_t> myStack;
//		std::stack<size_t> newOptimaStack;
//		std::vector<std::pair<size_t,size_t>> intervalsToClose;
//		myStack.push(nodeID);
//		while(!myStack.empty()){
//			size_t vertex=myStack.top();
//			bool descClosed=true;
//			if(candidateGraph.count(vertex)>0){
//				for(size_t desc:candidateGraph[vertex]){
//					if(isInRange(desc,boundLayer)){
//						if(isProcessed.count(desc)==0){
//							myStack.push(desc);
//							descClosed=false;
//						}
//					}
//				}
//
//				if(descClosed){
//					myStack.pop();
//					double bestValue=0;
//					std::unordered_set<size_t> bestNeighbors;
//					for(size_t desc:candidateGraph[vertex]){
//						//TODO also fill openForward!
//						if(valuesStructure[desc]<bestValue){ //TODO is contained in ValuesStr?
//							bestValue=valuesStructure;
//							bestNeighbors.clear();
//							bestNeighbors.insert(desc);
//						}
//						else if(valuesStructure[desc]==bestValue&&bestValue<0){
//							bestNeighbors.insert(desc);
//						}
//					}
//
//					if(isOneInOpt.count(vertex)>0){
//						for(size_t v:indexStructure[vertex]){
//							assert(bestNeighbors.count(v)>0);
//						}
//						for(size_t v:bestNeighbors){
//							if(indexStructure.count(v)==0){
//								newOptimaStack.push(v);
//								indexStructure[vertex].insert(v);
//								size_t obv=openBackward[vertex];
//								size_t ofv=openForward[v];
//								if(obv<ofv){
//									intervalsToClose.push_back({obv,ofv});
//								}
//								openBackward[v]=obv;
//								openForward.erase(v);
//								isOneInOpt.insert(v);
//
//							}
//						}
//					}
//					else{
//						size_t ofVertex=0;
//						for(size_t desc:bestNeighbors){
//							size_t descOf;
//							if(isOneInOpt.count(desc)>0){
//								descOf=openBackward[desc];
//							}
//							else{
//								descOf=openForward[desc];
//							}
//							if(descOf>ofVertex){
//								ofVertex=descOf;
//							}
//						}
//						openForward[vertex]=ofVertex;
//						indexStructure[vertex]=bestNeighbors;
//					}
//
//					//valuesStructure[vertex]=bestValue;
//					if(liftedCosts.count(vertex)>0){
//						bestValue+=liftedCosts[vertex];
//					}
//					valuesStructure[vertex]=bestValue; //TODO eventually delete bad values in valuesStructure and hence in candidateGraph
//				}
//			}
//
//		}
//
//		//Add all newly found optimal solutions
//		while(!newOptimaStack.empty()){
//			size_t newOptVertex=newOptimaStack.top();
//			newOptimaStack.pop();
//			isOneInOpt.insert(newOptVertex);
//			if(indexStructure.count(newOptVertex)>0){
//				for(size_t desc:indexStructure[newOptVertex]){
//					if(isOneInOpt.count(desc)==0){
//						openBackward[desc]=openBackward[newOptVertex];
//						openForward.erase(desc);
//						isOneInOpt.insert(desc);
//					}
//				}
//			}
//
//		}
//
//	}
//
//	//TODO close open vertices according to intervals to close
//	//TODO for all optimal vertices transform openBackward according to closed intervals
//	//TODO find candidates: for all optimal vertices, check their neighbors in candidateGraph (not optimal)
//	//TODO for found candidate, change open vertices in the identified interval, store messages
//
//
//
//	baseAndLiftedMessages messages={baseMessages,liftedMessages};
//	return messages;
//}





//template<class LDP_INSTANCE>
//inline baseAndLiftedMessages ldp_single_node_cut_factor<LDP_INSTANCE>::adjustCostsAndSendMessagesLifted(){
//	std::unordered_map<size_t,double> baseMessages=adjustCostsAndSendMessages();
//	std::unordered_map<size_t,double> liftedMessages;
//
//	std::unordered_set<size_t> isNotZeroInOpt;
//	std::unordered_set<size_t> isOneInOpt;
//	//TODO try other starting points then optimalSolution
//	std::stack<size_t> myStack;
//	myStack.push(optimalSolution);
//	std::vector<size_t> path;
//	std::unordered_set<size_t> isOnPath;
//	std::unordered_map<size_t,size_t> predOnPath;
//	bool pathClosed=false;
//
//
//	//TODO: Central vertex must be among the "parents of open vertices"
//
//	path.push_back(nodeID);
//	isOnPath.insert(nodeID);
//	size_t currentVertex=optimalSolution;
//	if(indexStructure[nodeID].size()>1){
//		for(size_t vertex:indexStructure[nodeID]){
//			if(vertex==optimalSolution)continue;
//			myStack.push(vertex);
//			predOnPath[vertex]=nodeID;
//		}
//	}
//	while(!pathClosed){
//		path.push_back(currentVertex);
//		isOnPath.insert(currentVertex);
//		isOneInOpt.insert(currentVertex);
//		isNotZeroInOpt.insert(currentVertex);
//		if(indexStructure.count(currentVertex)>0){
//			std::unordered_set<size_t>& myDesc=indexStructure[currentVertex];
//			if(myDesc.size()>1){
//				auto it=myDesc.begin();
//				it++;
//				for(;it!=myDesc.end();it++){
//					size_t sibling=*it;
//					myStack.push(sibling);
//					predOnPath[sibling]=currentVertex;
//				}
//			}
//			currentVertex=(*myDesc.begin());
//		}
//		else{
//			pathClosed=true;
//		}
//	}
//
//
//	size_t currentPathStart;  //does it need to be init?
//
//	std::unordered_map<size_t,std::unordered_set<size_t>> openNodesToParents;  //move above while, fill in in while
//	std::unordered_map<size_t,std::unordered_set<size_t>> parentsToConcurrentSiblings;
//
//	for (int i = 0; i < path.size()-1; ++i) {
//		openNodesToParents[path[i+1]].insert(path[i]);
//	}
//
//	while(!myStack.empty()){
//		size_t currentVertex=myStack.top();
//		myStack.pop();
//		if(predOnPath.count(currentVertex)>0){
//			currentPathStart=predOnPath[currentVertex];
//		}
//		isOneInOpt.insert(currentVertex);
//		bool hasDesc=indexStructure.count(currentVertex)>0;
//		if(!hasDesc){
//			size_t index=path.size()-1;
//
//			while(path[index]!=currentPathStart){
//				isNotZeroInOpt.erase(path[index]);
//				index--;
//			}
//		}
//		else{
//			std::unordered_set<size_t>& myDesc=indexStructure[currentVertex];
//			for(auto it=myDesc.begin();it!=myDesc.end();it++){
//				size_t vertex=*it;
//				if(isOnPath.count(vertex)>0){
//					size_t index=path.size()-1;
//					while(path[index]!=vertex){
//						index--;
//					}
//					index--;
//					while(path[index]!=currentPathStart){
//						isNotZeroInOpt.erase(path[index]);
//						if(openNodesToParents.count(path[index])>0){
//							for(size_t parent:openNodesToParents[path[index]]){
//								parentsToConcurrentSiblings.erase(parent);
//							}
//							openNodesToParents.erase(path[index]);
//						}
//						index--;
//					}
//					openNodesToParents[vertex].insert(currentVertex);
//					parentsToConcurrentSiblings[currentVertex]=std::unordered_set<size_t>();
//
//				}
//				else{
//					myStack.push(vertex);
//				}
//			}
//		}
//	}
//
//
//	std::unordered_map<size_t,std::unordered_set<size_t>> predecessors=initPredecessorsInIndexStr();
//	//Init the siblings here. They will disappear later. New should not be added to an existing parent later.
//	//Only new parents with new sets of siblings can be added later
//	for(auto it=openNodesToParents.begin();it!=openNodesToParents.end();it++){
//		size_t openVertex=it->first;
//		std::unordered_set<size_t>& parents=it->second;
//		double opt=valuesStructure[openVertex]; //Is it always in valuesStructure?
//		//if(-opt<minimalChange) minimalChange=-opt;
//		double concurrentValue=0;
//		//TODO: maybe here only create structure. Move number comparison to further steps
//		for(auto parent:parents){
//			for (int i = 0; i < numberOfNeighborsBase(parent); ++i) {
//				size_t sibling=getNeighborBaseVertex(parent,i);
//				if(sibling==openVertex) continue;
//				if(valuesStructure.count(sibling)>0){
//					double siblingValue=valuesStructure[sibling];
//					parentsToConcurrentSiblings[parent].insert(sibling);
//					//					if(siblingValue-opt<minimalChange){ //TODO if zero ->exception
//					//						minimalChange=siblingValue-opt;
//					//						minChangeNodeParSibling=std::vector<std::array<size_t,3>>(1);
//					//						minChangeNodeParSibling[0]={openVertex,parent,sibling};
//					//					}
//					//					else if(siblingValue-opt==minimalChange){
//					//						minChangeNodeParSibling.push_back({openVertex,parent,sibling});
//					//					}
//				}
//			}
//		}
//	}
//
//	while(!openNodesToParents.empty()){
//		double minimalChange=std::numeric_limits<double>::infinity();
//
//		std::unordered_map<size_t,std::array<size_t,2>> minChangeNodeParSibling;
//
//
//		for(auto it=openNodesToParents.begin();it!=openNodesToParents.end();it++){
//			size_t openVertex=it->first;
//			std::unordered_set<size_t>& parents=it->second;
//			double opt=valuesStructure[openVertex]; //Is it always in valuesStructure?
//			if(-opt<minimalChange){
//				minimalChange=-opt;
//				minChangeNodeParSibling[openVertex]={*(parents.begin()),getVertexToReach()};
//			}
//			double concurrentValue=0;   //TODO: default value corresponds to sibling==terminal vertex!
//			for(size_t parent:parents){
//				for(std::pair<size_t,double> siblingPair:parentsToConcurrentSiblings[parent]){  //maybe withou the double value, this will be updated
//					size_t sibling=siblingPair.first;
//					//if(sibling==openVertex) continue;  //This should not happen
//					if(valuesStructure.count(sibling)>0){  //This should be ensured!
//						double siblingValue=valuesStructure[sibling];  //should not be zero?
//						if(siblingValue-opt<minimalChange){ //TODO if zero ->exception
//							minimalChange=siblingValue-opt;
//							minChangeNodeParSibling.clear();
//							minChangeNodeParSibling[openVertex]={parent,sibling};
//						}
//						else if(siblingValue-opt==minimalChange){
//							minChangeNodeParSibling[openVertex]={parent,sibling};
//						}
//					}
//				}
//			}
//
//		}
//
//		size_t minChangeSize=minChangeNodeParSibling.size();
//		size_t vertexToClose;
//		if(minChangeSize>1){
//			//TODO select the more distant vertex
//			//vertexToClose=...
//		}
//		else{
//			vertexToClose=minChangeNodeParSibling.begin()->first;
//		}
//
//		std::stack<size_t> myStackNodes;
//		std::stack<double> myStackValues;
//		for(std::array<size_t,2> pair:minChangeNodeParSibling[vertexToClose]){
//			size_t sibling=pair[1];
//			myStackNodes.push(sibling);
//			if(liftedCosts.count(sibling)>0){
//				double cost=liftedCosts.at(sibling);
//				myStackValues.push(cost);
//			}
//			else{
//				myStackValues.push(0);
//			}
//		}
//		while(!myStack.empty()){
//			//DFS until reaching already optimal vertex
//			//Find all open vertices between vertex to close and the found vertex
//			//add respective value to every open node, store messages and close these nodes
//		}
//
//
//
//
//	}
//
//
//
//
//	//How to store parent and sibling information?
//
//
//
//
//
//
//
//}

//template<class LDP_INSTANCE>
//inline std::unordered_map<size_t,std::unordered_set<size_t>> ldp_single_node_cut_factor<LDP_INSTANCE>::initPredecessorsInIndexStr(){
//	std::unordered_map<size_t,std::unordered_set<size_t>> reachable;
//	std::unordered_map<size_t,std::unordered_set<size_t>> pred;
//	//TODO DFS from all vertices
//	for(auto pair:indexStructure){
//		std::size_t vertex=pair.first;
//		if(reachable.count(vertex)==0){
//			std::stack<size_t> stack;
//			stack.push(vertex);
//			while(!stack.empty()){
//				size_t currentVertex=stack.top();
//				if(indexStructure.count(currentVertex)==0){
//					reachable[currentVertex]=std::unordered_set<size_t>();
//					stack.pop();
//				}
//				else{
//					bool descClosed=true;
//					std::unordered_set<size_t>& myDesc=indexStructure[currentVertex];
//					for(size_t d:myDesc){
//						if(reachable.count(d)==0){
//							descClosed=false;
//							stack.push(d);
//						}
//					}
//					if(descClosed){
//						for(size_t d:myDesc){
//							reachable[currentVertex].insert(d);
//							reachable[currentVertex].insert(reachable[d].begin(),reachable[d].end()); //add all vertices reachable from desc.
//						}
//						for(size_t d:reachable[currentVertex]){
//							pred[d].insert(currentVertex);
//						}
//						stack.pop();
//					}
//				}
//			}
//		}
//	}
//	return pred;
//}

}
