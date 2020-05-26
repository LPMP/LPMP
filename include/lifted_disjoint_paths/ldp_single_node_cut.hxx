#pragma once
#include <andres/graph/digraph.hxx>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <array>
#include <list>
#include <set>
#include <config.hxx>

namespace LPMP {

//struct baseAndLiftedMessages{
//	std::unordered_map<size_t,double> baseMessages;
//	std::unordered_map<size_t,double> liftedMessages;
//};

struct StrForUpdateValues{
	std::unordered_map<size_t,double>& solutionCosts;
	std::unordered_map<size_t,double> valuesStructure;
	std::unordered_map<size_t,size_t> indexStructure;
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

	}
	bool useVertex(size_t vertex){
		if(useAllVertices) return true;
	    return valuesStructure.count(vertex)>0;
	}

	void setUseAllVertices(bool value){
		useAllVertices=value;
	}
	double getOptimalValue(){
		return optValue;
	}

};

template<class LDP_INSTANCE>
class ldp_single_node_cut_factor
{
public:
	constexpr static std::size_t nodeNotActive = std::numeric_limits<std::size_t>::max();

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


	bool hasLiftedEdgeToNode(size_t vertex) const{
		return liftedCosts.count(vertex)>0;
	}

	double LowerBound() const;

	double EvaluatePrimal() const;

	void init_primal(){
		primalBase_=nodeNotActive;
	}

	void setBaseEdgeActive(size_t vertex);
	void setNoBaseEdgeActive();

	void setPrimalLifted(std::unordered_set<size_t>& verticesOfActiveEdges);

	bool isNodeActive() const
	{
		return primalBase_!=nodeNotActive;
	}

	size_t getPrimalBase() const {
		return primalBase_;
	}

	const std::unordered_set<size_t>& getPrimalLifted() const {
		return primalLifted_;
	}

	const bool isActiveInPrimalLifted(size_t vertex) const {
		return primalLifted_.count(vertex)>0;
	}

	template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
	template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

	//auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }//TODO change this. This will not work with so many variables
	double tmp_to_delete_val;
	auto export_variables() { return std::tie(tmp_to_delete_val); } //?What comes here?

	void updateCostSimple(const double value,const size_t vertexIndex,bool isLifted);
	void updateNodeCost(const double value);
	double getOneBaseEdgeMinMarginal(size_t vertex);
	std::unordered_map<size_t,double> getAllBaseMinMarginals();

	std::unordered_map<size_t,double> getAllLiftedMinMarginals();

	double getNodeMinMarginal()const;
	double oneLiftedMinMarginal(size_t vertexOfLiftedEdge)const;



	const std::size_t nodeID;
	//const size_t nodeNotActive;
	size_t primalBase_;
	std::unordered_set<size_t> primalLifted_;





private:
	void updateValues() const;
	void updateValues(StrForUpdateValues& myStr,size_t vertexToIgnore=std::numeric_limits<size_t>::max()) const;
	std::unordered_map<size_t,double> bottomUpUpdate(const StrForUpdateValues& myStr,size_t vertex,std::unordered_map<size_t,size_t>& indexStr,std::unordered_set<size_t>* pClosedVert=0,std::unordered_map<size_t,double>* pBUValuesStr=0) const;
	void updateOptimal() const;

	std::list<size_t> getOptLiftedFromIndexStr(const StrForUpdateValues& myStr)const;


	std::unordered_map<size_t,std::unordered_set<size_t>> createCandidateGraph(const StrForUpdateValues& myStr);

//	struct vertexCompare {
//		vertexCompare(ldp_single_node_cut_factor<LDP_INSTANCE>& _sncFactor):sncFactor(_sncFactor){}
//		const ldp_single_node_cut_factor<LDP_INSTANCE>& sncFactor;
//	    bool operator() (const size_t& a, const size_t& b) const {
//	       return sncFactor.reachable(a,b);
//	    }
//	};
//	std::list<size_t>::iterator findAllOptimal(std::list<size_t>& isNotZeroInOpt,std::unordered_set<size_t>& isOneInOpt,std::unordered_map<size_t,std::unordered_set<size_t>>& candidateGraph,const StrForUpdateValues& strForUpdateValues);

	//std::set<size_t, decltype(vertexCompare)> s(vertexCompare);

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

	bool reachable(size_t firstVertex,size_t secondVertex)const{
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
	//mutable std::unordered_set<size_t> optimalSolutionLifted;
	mutable std::list<size_t> optimalSolutionLifted;
	//std::set<size_t> optimalSolutionLifted;

	//std::set<size_t,decltype(vertexCompare)> mySet;

	std::size_t minLayer;
	std::size_t maxLayer;

	const bool isOutFlow;

	const andres::graph::Digraph<>& baseGraph;
	const andres::graph::Digraph<>& liftedGraph;
	const LDP_INSTANCE& ldpInstance;


	std::unordered_map<size_t,double> baseCosts;
	std::unordered_map<size_t,double> liftedCosts;
    double nodeCost;
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
strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID)
//nodeNotActive(nID),strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID)
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
//	baseCosts[nodeNotActive]=0;
	optLiftedUpToDate=false;
	optBaseUpToDate=false;

	nodeCost=0;
	optValue=0;

}


template<class LDP_INSTANCE>
inline std::list<size_t> ldp_single_node_cut_factor<LDP_INSTANCE>::getOptLiftedFromIndexStr(const StrForUpdateValues& myStr) const{
	size_t vertexInOptimalPath=myStr.indexStructure.at(nodeID);

	std::list<size_t> optLifted;
//	std::cout<<"opt lifted: "<<std::endl;
	double optValueComputed=0;
	if(vertexInOptimalPath!=nodeNotActive){
		optValueComputed=myStr.baseCosts.at(vertexInOptimalPath);
	}

	bool hasOptDescendant=vertexInOptimalPath!=nodeNotActive;
	if(hasOptDescendant) optValueComputed+=nodeCost;
	while(hasOptDescendant){
	//	std::cout<<vertexInOptimalPath<<","<<std::endl;

		if(myStr.liftedCosts.count(vertexInOptimalPath)>0){
		//	std::cout<<"is lifted "<<std::endl;
			optLifted.push_back(vertexInOptimalPath);
			optValueComputed+=myStr.liftedCosts.at(vertexInOptimalPath);
		}
		auto it=myStr.indexStructure.find(vertexInOptimalPath);
		hasOptDescendant=it!=myStr.indexStructure.end();
		if(hasOptDescendant){
			vertexInOptimalPath=it->second;
		}
	}
	assert(std::abs(optValueComputed - myStr.optValue) <= 1e-8);
//	std::cout<<"opt value "<<optValueComputed<<std::endl;
//	std::cout<<"opt value in class "<<myStr.optValue<<std::endl;


	return optLifted;

}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateOptimal() const{
	if(!optLiftedUpToDate){

		//std::cout<<"updating lifted str."<<std::endl;

		updateValues();
	}
	else if(!optBaseUpToDate){
		//std::cout<<"updating base str."<<std::endl;
		optimalSolutionBase=nodeNotActive;
		double minValue=0;
		for (auto it=solutionCosts.begin();it!=solutionCosts.end();it++) {
			double value=it->second;
			if(value<minValue){
				minValue=value;
				optimalSolutionBase=it->first;
			}
		}
		strForUpdateValues.indexStructure[nodeID]=optimalSolutionBase;
		strForUpdateValues.optValue=minValue;
		optimalSolutionLifted=getOptLiftedFromIndexStr(strForUpdateValues);
		optBaseUpToDate=true;
	}
	optValue=strForUpdateValues.optValue;
}



template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setPrimalLifted(std::unordered_set<size_t>& verticesOfActiveEdges) {
	primalLifted_=verticesOfActiveEdges;
}




template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::EvaluatePrimal() const{
	//double value=0;
	if(primalBase_==nodeNotActive){
		return 0;
	}
	else{
		double value=nodeCost;
		value+=baseCosts.at(primalBase_);
		for(size_t node:primalLifted_){
			value+=liftedCosts.at(node);
		}
		return value;
	}
}




template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setBaseEdgeActive(size_t vertex){
	assert(vertex!=nodeID&&baseCosts.count(vertex)>0);
	primalBase_=vertex;

}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setNoBaseEdgeActive(){

	primalBase_=nodeNotActive;
}


template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getNodeMinMarginal()const{
	//std::cout<<"node min marginal "<<nodeID<<" "<<isOutFlow<<std::endl;
	updateOptimal();
	if(optimalSolutionBase!=nodeNotActive){
		return optValue-solutionCosts.at(nodeNotActive);
	}
	else{
		double value=std::numeric_limits<double>::infinity();
		for(auto pair:solutionCosts){
			if(pair.first==nodeNotActive) continue;
			if(pair.second<value){
				value=pair.second;
			}
		}
		return value-solutionCosts.at(nodeNotActive);
	}
}

template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getOneBaseEdgeMinMarginal(size_t vertex){
	assert(vertex!=nodeNotActive&&baseCosts.count(vertex)>0);
	updateOptimal();
	if(optimalSolutionBase!=vertex){
		return solutionCosts[vertex]-solutionCosts[optimalSolutionBase];
	}
	else{
		double secondBest=std::numeric_limits<double>::max();
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
	//std::cout<<"output min marginals"<<std::endl;
	std::unordered_map<size_t,double> minMarginals;
	if(optimalSolutionBase==nodeNotActive){
		double value=solutionCosts.at(nodeNotActive);
		for(auto pair:solutionCosts){
			if(pair.first!=nodeNotActive){
				minMarginals[pair.first]=pair.second-value;
				//std::cout<<pair.first<<": "<<(pair.second-value)<<std::endl;

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
				//std::cout<<it->first<<": "<<(value-secondBest)<<std::endl;
			}
		}
	}
	return minMarginals;
}

template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateNodeCost(const double value){
	//std::cout<<"Update node cost node: "<<nodeID<<" value "<<value<<std::endl;
	//double lbBefore=LowerBound();
	//std::cout<<"lb before"<<lbBefore<<std::endl;
	nodeCost+=value;
	//baseCosts[nodeNotActive]-=value;
	//solutionCosts[nodeNotActive]-=value;
	for(std::unordered_map<size_t,double>::iterator it=solutionCosts.begin();it!=solutionCosts.end();it++){
		if(it->first!=nodeNotActive){
		    it->second+=value;
		}
	}
	//double lbAfter=LowerBound();
//	std::cout<<"lb after"<<lbAfter<<std::endl;
	optBaseUpToDate=false;
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
		bool isOpt=false;
		for(size_t v:optimalSolutionLifted){
			if(v==vertexIndex){
				isOpt=true;
				break;
			}
		}
	//	std::cout<<"Update lifted cost "<<nodeID<<" "<<vertexIndex<<" value "<<value<<"optimal "<<isOpt<<std::endl;
		//double lbBefore=LowerBound();
		//std::cout<<"lb before"<<lbBefore<<std::endl;
		assert(liftedCosts.count(vertexIndex)>0);
		liftedCosts[vertexIndex]+=value;
		//valuesStructure[vertexIndex]+=value;
		optLiftedUpToDate=false;
		optBaseUpToDate=false;

	//	double lbAfter=LowerBound();
		//std::cout<<"lb after"<<lbAfter<<std::endl;
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
	//std::cout<<"update values run"<<std::endl;

	//StrForUpdateValues strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID);
	strForUpdateValues.indexStructure.clear();
	strForUpdateValues.solutionCosts.clear();
	strForUpdateValues.valuesStructure.clear();
	strForUpdateValues.setUseAllVertices(true);

	updateValues(strForUpdateValues);

	//std::cout<<"values updated"<<std::endl;
	optimalSolutionBase=strForUpdateValues.indexStructure[nodeID];
	optValue=solutionCosts[optimalSolutionBase];
	//std::cout<<"opt base: "<<optimalSolutionBase<<std::endl;

	optimalSolutionLifted=getOptLiftedFromIndexStr(strForUpdateValues);
	//std::cout<<std::endl;

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
		if(closedVertices.count(currentNode)>0){
			nodeStack.pop();
		}
		else{

			//std::cout<<"current vertex "<<currentNode<<std::endl;
			//std::cout<<"current node "<<currentNode<<std::endl;
			bool descClosed=true;
			double minValue=0;
			//std::unordered_set<size_t> minValueIndices;
			size_t minValueIndex=getVertexToReach();

			//std::cout<<"descendants: ";
			for (int i = 0; i < numberOfNeighborsBase(currentNode); ++i) {
				size_t desc=getNeighborBaseVertex(currentNode,i);
				if(desc==vertexToIgnore||desc==getVertexToReach()) continue;
				//std::cout<<desc;
				if(isInThisFactorRange(desc)&&myStr.useVertex(desc)){
					if(closedVertices.count(desc)>0||(lastLayerSet&&!isInGivenRange(desc,lastLayer))){  //descendant closed
						//std::cout<<"-cl";
						if(descClosed){
							auto it=myStr.valuesStructure.find(desc);
							if(it!=myStr.valuesStructure.end()){
								double value=it->second;
								if(minValue>value){
									minValue=value;
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
//				else{
//					std::cout<<"-nu";
//				}
//				std::cout<<", ";

			}
			//std::cout<<std::endl;
			if(descClosed){ //Close node if all descendants are closed
				//	std::cout<<"desc closed"<<std::endl;
				if(currentNode==nodeID){  //all nodes closed, compute solution values
					double bestValue=0;
					size_t bestVertex=nodeNotActive;
					myStr.solutionCosts[nodeNotActive]=0;
					for (auto it=myStr.baseCosts.begin();it!=myStr.baseCosts.end();++it) {
						//	std::cout<<"base edge "<<it->first<<std::endl;
						size_t vertex=it->first;
						if(vertex==nodeNotActive||vertex==vertexToIgnore) continue;
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
						double value=baseCost+nodeCost+valueToAdd;
						//std::cout<<"vertex "<<vertex<<", value: "<<value<<std::endl;
						myStr.solutionCosts[vertex]=value;
						if(value<bestValue){
							bestValue=value;
							bestVertex=vertex;
							//std::cout<<"best vertex"<<std::endl;
						}
						//std::cout<<"end for"<<std::endl;liftedCosts
					}
					//std::cout<<"after for"<<std::endl;
					//optimalSolution=bestVertex;
					//std::cout<<"final best vertex "<<bestVertex<<std::endl;
					myStr.indexStructure[nodeID]=bestVertex;
				//	std::cout<<nodeID<<"->"<<bestVertex<<std::endl;
					myStr.optValue=bestValue;
				}
				else{

					double valueToStore=minValue;
					assert(valueToStore<eps);
					const auto liftedIt=myStr.liftedCosts.find(currentNode);
					if(liftedIt!=myStr.liftedCosts.end()){
						valueToStore+=liftedIt->second;  //add lifted edge cost if the edge exists
					}
					const auto baseIt=baseCosts.find(currentNode);
					myStr.indexStructure[currentNode]=minValueIndex;
					//std::cout<<currentNode<<"->"<<minValueIndex<<", value: "<<valueToStore<<std::endl;

					if(valueToStore<0||(myStr.valuesStructure.count(currentNode)>0)||(baseIt!=myStr.baseCosts.end()&&minValue<0)){  //store only negative values or values needed to correct solutionCosts
						myStr.valuesStructure[currentNode]=valueToStore;
						//myStr.indexStructure[currentNode]=minValueIndex;
						//std::cout<<"vs insert"<<minValueIndex<<std::endl;
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


//template<class LDP_INSTANCE>
//inline std::list<size_t>::iterator ldp_single_node_cut_factor<LDP_INSTANCE>::findAllOptimal(std::list<size_t>& isNotZeroInOpt,std::unordered_set<size_t>& isOneInOpt,std::unordered_map<size_t,std::unordered_set<size_t>>& candidateGraph,const StrForUpdateValues& strForUpdateValues){
//	std::stack<size_t> myStack;
//	std::list<size_t> parentStack;
//
//	//TODO fix this function, return iterator to the beginning of the remaining list isNotZeroInOpt
//	std::unordered_set<size_t> closed;
//	myStack.push(nodeID);
//	parentStack.push_back(nodeID);
//	std::cout<<"find all optimal "<<std::endl;
//
//	//std::unordered_set<size_t> descendants;
//	double bestValue=strForUpdateValues.optValue;
//	for(auto pair:strForUpdateValues.solutionCosts){
//		size_t desc=pair.first;
//		double value=pair.second;
//		if(value==bestValue){
//			//descendants.insert(desc);
//			if(desc!=nodeNotActive&&desc!=getVertexToReach()){
//				myStack.push(desc);
//				closed.insert(desc);
//				std::cout<<"first in stack "<<desc<<std::endl;
//			}
//			else{
//				isNotZeroInOpt.clear();//TODO find out if this is necessary
//			}
//		}
//	}
////	for(size_t desc:strForUpdateValues.solutionCosts){
////		myStack.push(desc);
////	}
//
//	while(!myStack.empty()){
//		isOneInOpt.insert(myStack.top());
//		std::cout<<"processing vertex "<<myStack.top();
//		if(*(parentStack.rbegin())==myStack.top()){
//			std::cout<<"parent "<<myStack.top()<<" solved."<<std::endl;
//			parentStack.pop_back();
//			myStack.pop();
//
//		}
//		else{
//			size_t currentVertex=myStack.top();
//
//			double bestValue=0;
//			const auto descIt=strForUpdateValues.indexStructure.find(currentVertex);
//			if(descIt!=strForUpdateValues.indexStructure.end()){
//				size_t desc=descIt->second;
//				if(desc!=getVertexToReach()){
//					bestValue=strForUpdateValues.valuesStructure.at(desc);
//					std::cout<<"official best desc "<<desc<<std::endl;
//					if(closed.count(desc)==0){
//						myStack.push(desc);
//						closed.insert(desc);
//					}
//				}
//			}
//
//			std::unordered_set<size_t>& candidates=candidateGraph[currentVertex];
//			if(!candidates.empty()){
//				parentStack.push_back(currentVertex);
//				for(size_t desc:candidates){
//					if(desc==descIt->second) continue;
//
//					double value=strForUpdateValues.valuesStructure.at(desc);
//					if(value==bestValue){
//						std::cout<<"further opt desc "<<desc<<std::endl;
//						if(closed.count(desc)==0){
//							myStack.push(desc);
//							closed.insert(desc);
//						}
//					}
//
//				}
//			}
//			if(bestValue==0){  //The optimal path can be terminated in the current vertex
//				if(!isNotZeroInOpt.empty()){
//					std::unordered_set<size_t> keepOpen;
//					keepOpen.insert(currentVertex);
//					for(size_t optVertex:parentStack){
//						if(isNotZeroInOpt.count(optVertex)>0){
//							keepOpen.insert(optVertex);
//						}
//					}
//					isNotZeroInOpt=keepOpen;
//				}
//
//				myStack.pop();
//			}
//		}
//	}
//}



template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::oneLiftedMinMarginal(size_t vertexOfLiftedEdge)const{
	updateOptimal();


	bool hasEdge=liftedCosts.count(vertexOfLiftedEdge);
//	std::list<size_t> isNotZeroInOpt(optimalSolutionLifted.begin(),optimalSolutionLifted.end());
//	std::unordered_map<size_t,std::unordered_set<size_t>> candidateGraph=createCandidateGraph(strForUpdateValues);
//	findAllOptimal(isNotZeroInOpt,isOneInOpt,candidateGraph,strForUpdateValues);


	assert(hasEdge);
	bool isOptimal=false;
	for(size_t optVertex:optimalSolutionLifted){
		if(optVertex==vertexOfLiftedEdge){
			isOptimal=true;
			break;
		}
	}



	//std::cout<<"lilfted marginal "<<nodeID<<" "<<vertexOfLiftedEdge<<std::endl;

	if(isOptimal){
		//double lbBefore=LowerBound();
		//std::cout<<"OPTIMAL"<<std::endl;

		std::unordered_map<size_t,double> localSolutionCosts;
		StrForUpdateValues myStr(baseCosts,liftedCosts,localSolutionCosts,nodeID);
		myStr.indexStructure=strForUpdateValues.indexStructure;
		myStr.valuesStructure=strForUpdateValues.valuesStructure;
		myStr.optValue=strForUpdateValues.optValue;
		myStr.solutionCosts=strForUpdateValues.solutionCosts;


		myStr.setUseAllVertices(false);

		updateValues(myStr,vertexOfLiftedEdge);
		//size_t restrictedOptimalSolution=myStr.indexStructure[nodeID];
		double restrictedOptValue=myStr.optValue;

		double valueToReturn=optValue-restrictedOptValue;
		if(valueToReturn>1e-8){
			std::cout<<"wrong min marginal"<<nodeID<<", "<<vertexOfLiftedEdge<<", opt vertex but positive value "<<valueToReturn<<std::endl;
		}


		//std::cout<<"marginal "<<valueToReturn<<"node cost "<<nodeCost<<std::endl;
		return valueToReturn;


//		std::unordered_map<size_t,double> localSolutionCosts;
//		StrForUpdateValues myStr(strForUpdateValues);
//		myStr.useAllVertices=false;
//		updateValues(myStr,vertexOfLiftedEdge);
//		//size_t restrictedOptimalSolution=myStr.indexStructure[nodeID];
//		double restrictedOptValue=myStr.optValue;
//
//		return optValue-restrictedOptValue;

	}
	else{

		std::unordered_map<size_t,double> localLiftedCosts=liftedCosts;
		std::unordered_map<size_t,double> localSolCosts;
		std::unordered_map<size_t,double> localBaseCosts=baseCosts;
		StrForUpdateValues myStr(localBaseCosts,localLiftedCosts,localSolCosts,nodeID);

		updateValues(myStr);

		//TODO:ideally just use global structure strForUpdateValues in bottomUpUpdate
		std::unordered_map<size_t,size_t> indexStr;
		std::unordered_map<size_t,double> message=bottomUpUpdate(myStr,vertexOfLiftedEdge,indexStr);

		auto it =message.begin();
		double messValue=it->second;

//		StrForUpdateValues myStr(baseCosts,liftedCosts,solutionCosts,nodeID);

	//	StrForUpdateValues& myStr=strForUpdateValues;

//		strForUpdateValues.indexStructure.clear();
//		strForUpdateValues.solutionCosts.clear();
//		strForUpdateValues.valuesStructure.clear();
//		strForUpdateValues.setUseAllVertices(true);


	//	updateValues(strForUpdateValues);
	//	StrForUpdateValues tempStr=strForUpdateValues;

//		updateValues(myStr);

//		if(debug()){
//			for(auto pair:localLiftedCosts){
//				size_t vertex=pair.first;
//				double value=pair.second;
//				double value2=liftedCosts.at(vertex);
//				if(std::abs(value-value2)>eps){
//					std::cout<<"lifted cost mismatch "<<vertex<<": "<<value<<", "<<value2<<std::endl;
//				}
//			}
//		}
//
//		if(debug()){
//			size_t v=myStr.indexStructure[nodeID];
//			double optOrgCopmuted=0;
//
//			if(v!=nodeNotActive){
//				std::vector<size_t> solVertices;
//
//				optOrgCopmuted+=baseCosts.at(v);
//				optOrgCopmuted+=nodeCost;
//
//				do{
//					//std::cout<<v;
//					solVertices.push_back(v);
//					bool le=hasLiftedEdgeToNode(v);
//					if(le){
//						//	std::cout<<"(l)";
//						double cost=localLiftedCosts.at(v);
//						//std::cout<<cost;
//						optOrgCopmuted+=localLiftedCosts.at(v);
//					}
//					//std::cout<<", ";
//					auto it=myStr.indexStructure.find(v);
//					if(it!=myStr.indexStructure.end()){
//						v=it->second;
//					}
//					else break;
//				}while(v!=getVertexToReach());
//				//std::cout<<"opt orig computed "<<optOrgCopmuted<<std::endl;
//				if(v!=getVertexToReach()){
//					std::cout<<"last vertex is not vertex to reach "<<std::endl;
//				}
//
//
//				if(std::abs(optOrgCopmuted-myStr.optValue)>eps){
//					std::cout<<"mismatch between opt solutions "<<myStr.optValue<<", "<<optOrgCopmuted<<std::endl;
//					double cumulativeValue=0;
//					for(auto it=solVertices.rbegin();it!=solVertices.rend();it++){
//						size_t vert=*it;
//						std::cout<<vert;
//						auto liftIt=localLiftedCosts.find(vert);
//						if(liftIt!=localLiftedCosts.end()){
//							double lValue=liftIt->second;
//							cumulativeValue+=lValue;
//							std::cout<<"("<<lValue<<")";
//						}
//						std::cout<<std::endl;
//						double vsValue=myStr.valuesStructure[vert];
//						if(vsValue!=cumulativeValue){
//							std::cout<<"mismatch in cumulative values:  "<<cumulativeValue<<" "<<vsValue<<std::endl;
//						}
//					}
//					size_t firstVertex=solVertices[0];
//					cumulativeValue+=baseCosts.at(firstVertex);
//					cumulativeValue+=nodeCost;
//					double sc=myStr.solutionCosts.at(firstVertex);
//					if(cumulativeValue!=sc){
//						std::cout<<"mismatch in solution cost "<<cumulativeValue<<", "<<sc<<std::endl;
//					}
//				}
//			}
//		}






		//updateValues(myStr);
		//double origOptimal=myStr.optValue;
//		double origOptimal=myStr.optValue;
//		std::unordered_map<size_t,size_t> indexStr;
//		std::unordered_map<size_t,double> message=bottomUpUpdate(myStr,vertexOfLiftedEdge,indexStr);
//
//		auto it =message.begin();
//        double messValue=it->second;



//		localLiftedCosts[vertexOfLiftedEdge]-=messValue;
//
//
//		updateValues(myStr);
//
//		double newOptimal=myStr.optValue;
//		double eps=1e-8;
//		if(LPMP::debug()){
//			if(newOptimal+eps<origOptimal){
//
//				size_t v=strForUpdateValues.indexStructure[nodeID];
//				double optOrgCopmuted=0;
//				std::cout<<"recompute orig opt value"<<std::endl;
//				if(v!=nodeNotActive){
//					optOrgCopmuted+=baseCosts.at(v);
//					optOrgCopmuted+=nodeCost;
//				}
//				do{
//					std::cout<<v;
//					bool le=hasLiftedEdgeToNode(v);
//					if(le){
//						std::cout<<"(l)";
//						double cost=localLiftedCosts.at(v);
//						std::cout<<cost;
//						optOrgCopmuted+=localLiftedCosts.at(v);
//					}
//					std::cout<<", ";
//					auto it=strForUpdateValues.indexStructure.find(v);
//					if(it!=strForUpdateValues.indexStructure.end()){
//						v=it->second;
//					}
//					else break;
//				}while(v!=getVertexToReach());
//				std::cout<<"opt orig computed "<<optOrgCopmuted<<std::endl;
//				if(v!=getVertexToReach()){
//					std::cout<<"last vertex is not vertex to reach "<<std::endl;
//				}
//
//
//
//
//
//
//				double assumedOptValue=nodeCost;
//				double actualOptValue=nodeCost;
//
//				double newEdgeCost=localLiftedCosts[vertexOfLiftedEdge];
//				double oldLiftedEdgeCost=liftedCosts.at(vertexOfLiftedEdge);
			//	std::cout<<"wrong min marginal: "<<nodeID<<", "<<vertexOfLiftedEdge<<", not opt vertex, too high value of message: "<<(messValue)<<std::endl;
//				std::cout<<"new optimal "<<newOptimal<<", old optimal "<<origOptimal<<", new l. edge cost "<<newEdgeCost<<", old edge cost: "<<oldLiftedEdgeCost<<std::endl;
//				v=vertexOfLiftedEdge;
//				std::cout<<"opt sol. vertices from bottom up :";
//				do{
//					std::cout<<v<<", ";
//					auto liftIt=liftedCosts.find(v);
//					if(liftIt!=liftedCosts.end()){
//						assumedOptValue+=liftIt->second;
//					}
//					auto it=indexStr.find(v);
//
//					if(it!=indexStr.end()){
//						if(it->second==nodeID){
//							assumedOptValue+=baseCosts.at(v);
//						}
//						v=it->second;
//					}
//					else break;
//				}while(v!=nodeID);
//				if(v!=nodeID){
//					std::cout<<"last vertex is not "<<nodeID<<std::endl;
//				}
//
//
//				std::cout<<std::endl;
//				v=strForUpdateValues.indexStructure[vertexOfLiftedEdge];
//				do{
//					std::cout<<v<<", ";
//					auto liftIt=liftedCosts.find(v);
//					if(liftIt!=liftedCosts.end()){
//						assumedOptValue+=liftIt->second;
//					}
//					auto it=strForUpdateValues.indexStructure.find(v);
//					if(it!=strForUpdateValues.indexStructure.end()){
//						v=it->second;
//					}
//					else break;
//				}while(v!=getVertexToReach());
//				if(v!=getVertexToReach()){
//					std::cout<<"last vertex is not vertex to reach "<<std::endl;
//				}
//
//				assumedOptValue-=messValue;
//
//				std::cout<<std::endl;
//				std::cout<<"assumed opt value "<<assumedOptValue<<std::endl;
//				std::cout<<"opt from new opt run: ";
//				double optFromNewRun=0;
//				v=myStr.indexStructure[nodeID];
//
//				if(v!=nodeNotActive){
//					optFromNewRun+=baseCosts.at(v);
//					optFromNewRun+=nodeCost;
//				}
//				do{
//					std::cout<<v;
//					bool le=hasLiftedEdgeToNode(v);
//					if(le){
//						std::cout<<"(l)";
//						double cost=localLiftedCosts.at(v);
//						std::cout<<cost;
//						optFromNewRun+=localLiftedCosts.at(v);
//					}
//					std::cout<<", ";
//					auto it=myStr.indexStructure.find(v);
//					if(it!=myStr.indexStructure.end()){
//						v=it->second;
//					}
//					else break;
//				}while(v!=getVertexToReach());
//				std::cout<<"opt from new run computed "<<optFromNewRun<<std::endl;
//				if(v!=getVertexToReach()){
//					std::cout<<"last vertex is not vertex to reach "<<std::endl;
//				}
//				std::cout<<std::endl;
//				std::cout<<"original opt solution"<<std::endl;
//				for(size_t optVertex:optimalSolutionLifted){
//					std::cout<<optVertex<<", ";
//				}
//				std::cout<<std::endl<<std::endl;
//
//			}
//
//
//
//			if(messValue<-eps){
//				std::cout<<"wrong min marginal"<<nodeID<<", "<<vertexOfLiftedEdge<<", not opt vertex but negative value "<<messValue<<std::endl;
		//	}
	//	}
		return messValue;
	}

}


template<class LDP_INSTANCE>
inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::bottomUpUpdate(const StrForUpdateValues& myStr,size_t vertex,std::unordered_map<size_t,size_t>& indexStr,std::unordered_set<size_t>* pClosedVert,std::unordered_map<size_t,double>* pBUValuesStr)const{
	bool onlyOne=pClosedVert==0;

	std::unordered_map<size_t,double> messages;
	if(onlyOne){
		pClosedVert=new std::unordered_set<size_t> ();
		pBUValuesStr=new std::unordered_map<size_t,double>();
	}
	std::unordered_set<size_t>& closedVertices=*pClosedVert;
	std::unordered_map<size_t,double>& buValuesStr=*pBUValuesStr;

	std::stack<size_t> myStack;
	myStack.push(vertex);
	while(!myStack.empty()){
		size_t currentVertex=myStack.top();
		if(closedVertices.count(currentVertex)>0){
			myStack.pop();
		}
		else{
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
				size_t bestIndex=nodeID;
				double bestValue=std::numeric_limits<double>::infinity();
				auto baseIt=myStr.baseCosts.find(currentVertex);
				if(baseIt!=myStr.baseCosts.end()){
					bestValue=(baseIt->second)+nodeCost;
				}
				for (int i = 0; i < numberOfNeighborsBaseRev(currentVertex); ++i) {
					size_t pred=getNeighborBaseVertexRev(currentVertex,i);
					if(pred==nodeID) continue;
					auto valuesIt=buValuesStr.find(pred);
					if(valuesIt!=buValuesStr.end()){
						double value=valuesIt->second;
						if(value<bestValue){
							bestValue=value;
							bestIndex=valuesIt->first;
						}
					}
				}
				std::unordered_map<size_t,double>::const_iterator liftedIt=myStr.liftedCosts.find(currentVertex);
				indexStr[currentVertex]=bestIndex;
				if(liftedIt!=myStr.liftedCosts.end()){
					if(onlyOne&&currentVertex!=vertex){
						bestValue+=liftedIt->second;
					}
					else{
						bestValue+=liftedIt->second;
						double topDownValue=0;
						auto bestTdIt=myStr.indexStructure.find(currentVertex);
						size_t bestTd=bestTdIt->second;
//						if(liftedCosts.count(bestTd)>0&&isOneInOpt.count(bestTd)==0){
//
//							//std::cout<<"new optimal vertex "<<bestTd<<std::endl;
//
//						}
						if(bestTd!=getVertexToReach()){
							topDownValue=myStr.valuesStructure.at(bestTd);
						}

						//					auto valuesIt=myStr.valuesStructure.find(currentVertex);
						//					if(valuesIt!=myStr.valuesStructure.end()){
						//						topDownValue=valuesIt->second;
						//					}
						double restrictedOpt=topDownValue+bestValue;


						double delta=restrictedOpt-myStr.optValue;

						//std::cout<<"message "<<currentVertex<<": "<<delta<<std::endl;

						messages[currentVertex]=delta;
						bestValue-=delta;

//						size_t v=myStr.indexStructure.at(currentVertex);
//						double checkTDValue=0;
//						while(currentVertex!=getVertexToReach()){
//							auto lIt=myStr.liftedCosts.find(v);
//							if(lIt!=myStr.liftedCosts.end()){
//								checkTDValue+=lIt->second;
//							}
//							v=myStr.indexStructure.at(v);
//						}
//						double checkBUValue=0;
//						v=currentVertex;
//						while(currentVertex!=getVertexToReach()){
//							auto lIt=myStr.liftedCosts.find(v);
//							if(lIt!=myStr.liftedCosts.end()){
//								checkTDValue+=lIt->second;
//							}
//							v=myStr.indexStructure.at(v);
//						}
//

						//liftedIt->second-=delta;  //cannot be, lifted edge costs in myStr are const

					}
				}

				closedVertices.insert(currentVertex);
				buValuesStr[currentVertex]=bestValue;
				myStack.pop();

			}
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
	std::unordered_set<size_t> isClosed;
	for(auto it=myStr.valuesStructure.begin();it!=myStr.valuesStructure.end();++it){
		size_t vertex=it->first;
		if(liftedCosts.count(vertex)==0&&isClosed.count(vertex)==0){
			std::stack<size_t> nodeStack;
			nodeStack.push(vertex);
			while(!nodeStack.empty()){
				size_t currentVertex=nodeStack.top();
				bool descClosed=true;
				std::unordered_set<size_t> onlyLiftedDesc;
				for (int i = 0; i < numberOfNeighborsBase(currentVertex); ++i) {
					size_t neighbor=getNeighborBaseVertex(currentVertex,i);
					if(myStr.valuesStructure.count(neighbor)>0){
						if(liftedCosts.count(neighbor)==0){
							if(isClosed.count(neighbor)==0){
								descClosed=false;
								nodeStack.push(neighbor);

							}
							else if(descClosed){
								std::unordered_set<size_t> &neighborLiftedDesc=candidateGraph[neighbor];
								onlyLiftedDesc.insert(neighborLiftedDesc.begin(),neighborLiftedDesc.end());
							}
						}
						else if(descClosed){
							onlyLiftedDesc.insert(neighbor);
						}
					}
				}
				if(descClosed){
					candidateGraph[currentVertex]=onlyLiftedDesc;
					nodeStack.pop();
					isClosed.insert(currentVertex);
				}
			}
		}
	}

	std::unordered_map<size_t,std::unordered_set<size_t>> finalCandidateGraph;

	for(auto it=myStr.valuesStructure.begin();it!=myStr.valuesStructure.end();++it){
		size_t vertex=it->first;
		if(liftedCosts.count(vertex)>0){
			for (int i = 0; i < numberOfNeighborsBase(vertex); ++i) {
				size_t neighbor=getNeighborBaseVertex(vertex,i);
				if(myStr.valuesStructure.count(neighbor)>0){
					if(liftedCosts.count(neighbor)>0){
						finalCandidateGraph[vertex].insert(neighbor);
					}
					else{
						std::unordered_set<size_t>& neighborsDesc=candidateGraph[neighbor];
						finalCandidateGraph[vertex].insert(neighborsDesc.begin(),neighborsDesc.end());
					}
				}
			}
		}
	}

//		for(auto it=myStr.valuesStructure.begin();it!=myStr.valuesStructure.end();++it){
//			size_t vertex=it->first;
//			//strForUpdateValues.relevantVertices.insert(vertex);
//
//			for (int i = 0; i < numberOfNeighborsBase(vertex); ++i) {
//				size_t neighbor=getNeighborBaseVertex(vertex,i);
//				if(myStr.valuesStructure.count(neighbor)>0){
//					candidateGraph[vertex].insert(neighbor);
//				}
//			}
//		}

		return finalCandidateGraph;
}

template<class LDP_INSTANCE>
inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllLiftedMinMarginals(){


	updateOptimal();
	std::unordered_map<size_t,double> liftedMessages;

	double currentOptValue=strForUpdateValues.optValue;

	std::list<size_t> isNotZeroInOpt=optimalSolutionLifted;
	std::unordered_set<size_t> isOneInOpt(isNotZeroInOpt.begin(),isNotZeroInOpt.end());

	std::unordered_map<size_t,double> localSolutionCosts=solutionCosts;
	std::unordered_map<size_t,double> localLiftedCosts=liftedCosts;
	std::unordered_map<size_t,double> localBaseCosts=baseCosts;

	StrForUpdateValues myStr(localBaseCosts,localLiftedCosts,localSolutionCosts,nodeID);
	myStr.indexStructure=strForUpdateValues.indexStructure;
	myStr.valuesStructure=strForUpdateValues.valuesStructure;
	myStr.optValue=strForUpdateValues.optValue;
	myStr.solutionCosts=strForUpdateValues.solutionCosts;


	myStr.setUseAllVertices(false);


	auto listIt=isNotZeroInOpt.begin();
	while(!isNotZeroInOpt.empty()){
		size_t vertexToClose=*listIt;

		//std::cout<<"vertex to close "<<vertexToClose<<std::endl;
		updateValues(myStr,vertexToClose);
		double newOpt=myStr.optValue;
		//std::cout<<"new opt "<<newOpt<<std::endl;

		std::list<size_t> secondBest=getOptLiftedFromIndexStr(myStr);
		auto sbIt=secondBest.begin();

		listIt=isNotZeroInOpt.erase(listIt);
		while(listIt!=isNotZeroInOpt.end()&&sbIt!=secondBest.end()){
			if(*sbIt==*listIt){
				isOneInOpt.insert(*sbIt);
				sbIt++;
				listIt++;
			}
			else if(reachable(*sbIt,*listIt)){
				isOneInOpt.insert(*sbIt);
				sbIt++;
			}
			else if(reachable(*listIt,*sbIt)){
				listIt=isNotZeroInOpt.erase(listIt);
			}
			else{
				listIt=isNotZeroInOpt.erase(listIt);
				isOneInOpt.insert(*sbIt);
				sbIt++;
			}
		}
		isNotZeroInOpt.erase(listIt,isNotZeroInOpt.end());
		while(sbIt!=secondBest.end()){
			isOneInOpt.insert(*sbIt);
			sbIt++;
		}

		listIt=isNotZeroInOpt.begin();

		double delta=currentOptValue-newOpt;
		//std::cout<<"orig lifted cost "<<myStr.liftedCosts.at(vertexToClose)<<std::endl;
		localLiftedCosts[vertexToClose]-=delta;
		liftedMessages[vertexToClose]=delta;
		currentOptValue=newOpt;

		//std::cout<<"message "<<vertexToClose<<": "<<delta<<std::endl;
		//std::cout<<"delta for "<<vertexToClose<<": "<<delta<<", new l.cost: "<<localLiftedCosts[vertexToClose]<<std::endl;
		//std::cout<<"lifted cost in myStr "<<myStr.liftedCosts.at(vertexToClose)<<std::endl;
		//listIt=isNotZeroInOpt.erase(listIt);

	}

	myStr.setUseAllVertices(true);

	updateValues(myStr);
	std::unordered_map<size_t,double> buValuesStructure;
	std::unordered_set<size_t> closedVertices;
	for(size_t optVertex:isOneInOpt){
		buValuesStructure[optVertex]=currentOptValue-myStr.valuesStructure[optVertex]+localLiftedCosts[optVertex];
		closedVertices.insert(optVertex);
	}
	std::unordered_map<size_t,size_t> indexStr;
	for(auto pair:liftedCosts){
		if(closedVertices.count(pair.first)==0){

			std::unordered_map<size_t,double> newMessages=bottomUpUpdate(myStr,pair.first,indexStr,&closedVertices,&buValuesStructure);
			liftedMessages.insert(newMessages.begin(),newMessages.end());
		}
	}


	return liftedMessages;
}





class ldp_snc_lifted_message
{
public:
	ldp_snc_lifted_message(const std::size_t _left_node, const std::size_t _right_node)
	: left_node(_left_node),
	  right_node(_right_node)
	{}

	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamLeft(SINGLE_NODE_CUT_FACTOR& l, const double msg, const std::size_t msg_dim) const
	{
		assert(msg_dim == 0);
		l.updateCostSimple(msg,right_node,true);
	}

	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
	{
		assert(msg_dim == 0);
		r.updateCostSimple(msg,left_node,true);
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
	{
		const double delta = r.oneLiftedMinMarginal(left_node);
		msg[0] -= omega * delta;
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_right(const SINGLE_NODE_CUT_FACTOR& l, MSG& msg, const double omega)
	{
		const double delta = l.oneLiftedMinMarginal(right_node);
		msg[0] -= omega * delta;
	}

	template<typename SINGLE_NODE_CUT_FACTOR>
	bool check_primal_consistency(const SINGLE_NODE_CUT_FACTOR& l, const SINGLE_NODE_CUT_FACTOR& r) const
	{
		const bool left_snc_edge = l.isActiveInPrimalLifted(right_node);
		const bool right_snc_edge = r.isActiveInPrimalLifted(left_node);
		return left_snc_edge == right_snc_edge;
	}

private:
	std::size_t left_node;
	std::size_t right_node;
};


class ldp_snc_node_message
{
public:
	ldp_snc_node_message(const std::size_t _node)
	: node(_node)
	{}

	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamLeft(SINGLE_NODE_CUT_FACTOR& l, const double msg, const std::size_t msg_dim) const
	{
		assert(msg_dim == 0);
		l.updateNodeCost(msg);
	}

	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
	{
		assert(msg_dim == 0);
		r.updateNodeCost(msg);
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
	{
		const double delta = r.getNodeMinMarginal();
		msg[0] -= omega * delta;
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_right(const SINGLE_NODE_CUT_FACTOR& l, MSG& msg, const double omega)
	{
		const double delta = l.getNodeMinMarginal();
		msg[0] -= omega * delta;
	}

	template<typename SINGLE_NODE_CUT_FACTOR>
	bool check_primal_consistency(const SINGLE_NODE_CUT_FACTOR& l, const SINGLE_NODE_CUT_FACTOR& r) const
	{
		const bool left_snc_active = l.isNodeActive();
		const bool right_snc_active = r.isNodeActive();
		return left_snc_active == right_snc_active;
	}

private:
	std::size_t node;
};

}
