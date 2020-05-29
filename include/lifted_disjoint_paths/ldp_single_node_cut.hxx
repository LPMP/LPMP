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
	std::vector<double>& solutionCosts;
	std::unordered_map<size_t,double> valuesStructure;
	std::unordered_map<size_t,size_t> vertexIDStructure;
	size_t optBase;
	const std::vector<double>& baseCosts;
	const std::vector<double>& liftedCosts;



	//std::unordered_set<size_t> relevantVertices;
	bool useAllVertices;
	double optValue;
	const size_t nodeID;
	//size_t optimalSolution;
	StrForUpdateValues(const std::vector<double>& bCosts,const std::vector<double>& lCosts,std::vector<double>& solCosts,const size_t centralNode):
	baseCosts(bCosts),
	liftedCosts(lCosts),
	solutionCosts(solCosts),
	optValue(0),
	nodeID(centralNode)
	{
		useAllVertices=true;
		optBase=0;

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

	void copyFromOther(StrForUpdateValues& str){
		vertexIDStructure=str.vertexIDStructure;
		valuesStructure=str.valuesStructure;
		solutionCosts=str.solutionCosts;
		optValue=str.optValue;
	}

};

template<class LDP_INSTANCE>
class ldp_single_node_cut_factor
{
public:
	//constexpr static std::size_t nodeNotActive = std::numeric_limits<std::size_t>::max();

	//template<class LPD_STRUCT> ldp_single_node_cut_factor(const LPD_STRUCT& ldpStruct);
	ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut);
	ldp_single_node_cut_factor(const ldp_single_node_cut_factor&);

	void initBaseCosts(double fractionBase);

	void initLiftedCosts(double fractionLifted);

	//void addLiftedEdge(size_t node,double cost);

	const andres::graph::Digraph<>& getBaseGraph() const {
		return baseGraph;
	}

	const std::vector<double>& getLiftedCosts() const {
		return liftedCosts;
	}

	const std::vector<double>& getBaseCosts() const {
		return baseCosts;
	}


//	bool hasLiftedEdgeToNode(size_t vertex) const{
//		return liftedCosts.count(vertex)>0;
//	}

	double LowerBound() const;

	double EvaluatePrimal() const;

	void init_primal(){
		primalBase_=nodeNotActive;
	}

	void setBaseEdgeActive(size_t index);
	void setNoBaseEdgeActive();

	void setPrimalLifted(std::unordered_set<size_t>& verticesOfActiveEdges);

	bool isNodeActive() const
	{
		return primalBase_!=nodeNotActive;
	}

	size_t getPrimalBaseIndex() const {
		return primalBase_;
	}

	size_t getPrimalBaseVertexID() const {
		return baseIDs[primalBase_];
	}

	const std::unordered_set<size_t>& getPrimalLiftedIndices() const {
		return primalLifted_;
	}

	const bool isActiveInPrimalLifted(size_t vertex) const {
		return primalLifted_.count(vertex)>0;
		//return primalLifted_.at(vertex);
	}

	const size_t getLiftedOrderToID(size_t vertexID) const{
		return liftedIDToOrder.at(vertexID);
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

	std::vector<double> getAllLiftedMinMarginals();

	double getNodeMinMarginal()const;
	double oneLiftedMinMarginal(size_t indexOfLiftedEdge)const;



	const std::size_t nodeID;
	//const size_t nodeNotActive;
	size_t primalBase_;
	std::unordered_set<size_t> primalLifted_;

	//TODO use this
	//std::vector<bool> primalLifted_;

	std::vector<size_t> baseIDs;
    std::vector<size_t> liftedIDs;
	size_t nodeNotActive;




private:
	void updateValues() const;
	void updateValues(StrForUpdateValues& myStr,size_t vertexIDToIgnore) const;
	void updateValues(StrForUpdateValues& myStr) const;
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
	bool isInThisFactorRange(const size_t vertexID) const {
		assert(vertexID < baseGraph.numberOfVertices());
		if(isOutFlow){
			//if(debug()) std::cout<<"vertex: "<<vertexID<<", max l. "<<maxLayer<<", layer of vertex "<<ldpInstance.getGroupIndex(vertexID);
			if(vertexID==ldpInstance.getTerminalNode()) return true;
			else return ldpInstance.getGroupIndex(vertexID)<=maxLayer;
		}
		else{
			if(vertexID==ldpInstance.getSourceNode()) return true;
			else return ldpInstance.getGroupIndex(vertexID)>=minLayer;
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
	mutable std::list<size_t> optimalSolutionLifted;  //maybe change to vector<bool>

	//std::set<size_t> optimalSolutionLifted;

	//std::set<size_t,decltype(vertexCompare)> mySet;

	std::size_t minLayer;
	std::size_t maxLayer;

	const bool isOutFlow;

	const andres::graph::Digraph<>& baseGraph;
	const andres::graph::Digraph<>& liftedGraph;
	const LDP_INSTANCE& ldpInstance;
	size_t nodeNotActiveForStructures;


//	std::unordered_map<size_t,double> baseCosts;
//	std::unordered_map<size_t,double> liftedCosts;

	std::vector<double> baseCosts;
	std::vector<double> liftedCosts;
    double nodeCost;
	mutable std::vector<double> solutionCosts;



	std::unordered_map<size_t,size_t> baseIDToOrder;
	std::unordered_map<size_t,size_t> liftedIDToOrder;

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
inline  ldp_single_node_cut_factor<LDP_INSTANCE>::ldp_single_node_cut_factor(const ldp_single_node_cut_factor& sncFactor):
baseGraph(sncFactor.baseGraph),
liftedGraph(sncFactor.liftedGraph),
nodeID(sncFactor.nodeID),
ldpInstance(sncFactor.ldpInstance),
isOutFlow(sncFactor.isOutFlow),
baseCosts(sncFactor.baseCosts),
liftedCosts(sncFactor.liftedCosts),
solutionCosts(sncFactor.solutionCosts),
strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID),
nodeNotActive(sncFactor.nodeNotActive)
{


	primalBase_=sncFactor.primalBase_;
	optimalSolutionBase=sncFactor.optimalSolutionBase;
	minLayer=sncFactor.minLayer;
	maxLayer=sncFactor.maxLayer;

	nodeCost=sncFactor.nodeCost;
	optValue=sncFactor.optValue;

	strForUpdateValues.copyFromOther(sncFactor.strForUpdateValues);


	optLiftedUpToDate=sncFactor.optLiftedUpToDate;
	optBaseUpToDate=sncFactor.optBaseUpToDate;

	tmp_to_delete_val=sncFactor.tmp_to_delete_val;
	nodeNotActiveForStructures=sncFactor.nodeNotActiveForStructures;

}


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
	//if(debug()) std::cout<<" snc factor "<<nodeID<<", "<<isOut<<std::endl;


	if(isOutFlow){
		minLayer=ldpInst.getGroupIndex(nodeID);
		//if(debug()) std::cout<<"gap lifted "<<ldpInst.getGapLifted()<<std::endl;
		maxLayer=minLayer+ldpInst.getGapLifted(); //some method that returns max time gap lifted
		//if(debug()) std::cout<<"max layer "<<maxLayer<<std::endl;
	}
	else{
		maxLayer=ldpInst.getGroupIndex(nodeID);
		minLayer=std::max(0,int(maxLayer)-int(ldpInst.getGapLifted()));
	}

//	if(debug()) std::cout<<"min and max layers set "<<std::endl;
	//baseCosts=std::unordered_map<size_t,double>();
	initBaseCosts(0);
	//if(debug()) std::cout<<"base costs set "<<std::endl;
	//liftedCosts=std::unordered_map<size_t,double>();
	initLiftedCosts(0);
	//if(debug()) std::cout<<"lifted costs set "<<std::endl;
	nodeNotActive=baseCosts.size();
	//solutionCosts=std::unordered_map<size_t,double>();
	//solutionCosts[nodeNotActive]=0;
//	baseCosts[nodeNotActive]=0;
	optLiftedUpToDate=false;
	optBaseUpToDate=false;

	primalBase_=nodeNotActive;  //corresponds to no edge active
	optimalSolutionBase=nodeNotActive;
	nodeCost=0;
	optValue=0;
	nodeNotActiveForStructures=nodeID;




}


template<class LDP_INSTANCE>
inline std::list<size_t> ldp_single_node_cut_factor<LDP_INSTANCE>::getOptLiftedFromIndexStr(const StrForUpdateValues& myStr) const{

	std::list<size_t> optLifted;
   if(debug())	std::cout<<"opt lifted: "<<std::endl;
	double optValueComputed=0;
	if(myStr.optBase!=nodeNotActive){
		optValueComputed=myStr.baseCosts.at(myStr.optBase);
		optValueComputed+=nodeCost;

		bool hasOptDescendant=true;
		size_t vertexInOptimalPath=myStr.vertexIDStructure.at(baseIDs[myStr.optBase]);

		while(hasOptDescendant){
			if(debug())	std::cout<<vertexInOptimalPath<<","<<std::endl;

			if(liftedIDToOrder.count(vertexInOptimalPath)>0){
				//	std::cout<<"is lifted "<<std::endl;
				optLifted.push_back(vertexInOptimalPath);
				if(debug()) optValueComputed+=myStr.liftedCosts[liftedIDToOrder.at(vertexInOptimalPath)];
			}
			auto it=myStr.vertexIDStructure.find(vertexInOptimalPath);
			hasOptDescendant=it!=myStr.vertexIDStructure.end();
			if(hasOptDescendant){
				vertexInOptimalPath=it->second;
			}
		}
	}
	if(debug())	assert(std::abs(optValueComputed - myStr.optValue) <= 1e-8);
//	std::cout<<"opt value "<<optValueComputed<<std::endl;
//	std::cout<<"opt value in class "<<myStr.optValue<<std::endl;


	return optLifted;

}




//template<class LDP_INSTANCE>
//inline std::list<size_t> ldp_single_node_cut_factor<LDP_INSTANCE>::getOptLiftedFromIndexStr(const StrForUpdateValues& myStr) const{
//	size_t vertexInOptimalPath=myStr.indexStructure.at(nodeID);
//
//	std::list<size_t> optLifted;
////	std::cout<<"opt lifted: "<<std::endl;
//
//	std::unordered_set<size_t> optVertices;
//
//	bool hasOptDescendant=vertexInOptimalPath!=nodeNotActive;
//	while(hasOptDescendant){
//        optVertices.insert(vertexInOptimalPath);
//		auto it=myStr.indexStructure.find(vertexInOptimalPath);
//		hasOptDescendant=it!=myStr.indexStructure.end();
//		if(hasOptDescendant){
//			vertexInOptimalPath=it->second;
//		}
//	}
//
//	for (int i = 0; i < liftedIDs.size(); ++i) {
//		if(optVertices.count(liftedIDs[i])>0){
//
//		}
//	}
////	std::cout<<"opt value "<<optValueComputed<<std::endl;
////	std::cout<<"opt value in class "<<myStr.optValue<<std::endl;
//
//
//	return optLifted;
//
//}


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
		for (int i = 0; i < solutionCosts.size()-1; ++i) {
			double value=solutionCosts[i];
			if(value<minValue){
				minValue=value;
				optimalSolutionBase=i;
			}
		}


		if(optimalSolutionBase)
		strForUpdateValues.optBase=optimalSolutionBase;
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
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setBaseEdgeActive(size_t index){
	assert(index<baseCosts.size());
	primalBase_=index;

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
		for (int i = 0; i < solutionCosts.size()-1; ++i) {
			if(solutionCosts[i]<value){
				value=solutionCosts[i];
			}
		}

		return value-solutionCosts.at(nodeNotActive);
	}
}

template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getOneBaseEdgeMinMarginal(size_t index){
	assert(index<baseCosts.size());
	updateOptimal();
	if(optimalSolutionBase!=index){
		return solutionCosts[index]-solutionCosts[optimalSolutionBase];
	}
	else{
		double secondBest=std::numeric_limits<double>::max();
		for (int i = 0; i < solutionCosts.size(); ++i) {
			if(i==optimalSolutionBase) continue;
			if(solutionCosts[i]<secondBest){
				secondBest=solutionCosts[i];
			}
		}
		return solutionCosts[index]-secondBest;
	}
}

//TODO make this function return a vector
template<class LDP_INSTANCE>
inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllBaseMinMarginals(){
	updateOptimal();
	if(debug()) std::cout<<"output min marginals"<<std::endl;
	std::unordered_map<size_t,double> minMarginals;
	if(optimalSolutionBase==nodeNotActive){
		for (int i = 0; i < solutionCosts.size()-1; ++i) {
			minMarginals[i]=solutionCosts[i];
			if(debug()) std::cout<<i<<" "<<minMarginals[i]<<std::endl;

		}
	}
	else{
		double secondBest=std::numeric_limits<double>::infinity();
		double optValue=solutionCosts.at(optimalSolutionBase);
		for (int i = 0; i < solutionCosts.size(); ++i) {
			if(i==optimalSolutionBase) continue;
			if(solutionCosts[i]<secondBest){
				secondBest=solutionCosts[i];
			}

		}
		for (int i = 0; i < solutionCosts.size()-1; ++i) {
			minMarginals[i]=solutionCosts[i]-secondBest;
			if(debug()) std::cout<<i<<" "<<minMarginals[i]<<std::endl;
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
	for (int i = 0; i < solutionCosts.size()-1; ++i) {
		(solutionCosts[i]+=value);
	}

	optBaseUpToDate=false;
}

template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateCostSimple(const double value,const size_t vertexIndex,bool isLifted){//Only cost change
	if(!isLifted){ //update in base edge
		assert(vertexIndex<baseCosts.size());
		baseCosts[vertexIndex]+=value;
		solutionCosts[vertexIndex]+=value;
		optBaseUpToDate=false;
	}
	else{ //update in lifted edge

		assert(vertexIndex<liftedCosts.size());
		liftedCosts[vertexIndex]+=value;

		optLiftedUpToDate=false;
		optBaseUpToDate=false;


	}
}


//template<class LDP_INSTANCE>
//inline void ldp_single_node_cut_factor<LDP_INSTANCE>::addLiftedEdge(size_t node,double cost){
//	assert(reachable(nodeID,node));
//	liftedCosts[node]=cost;
//}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::initLiftedCosts(double fractionLifted){
	liftedCosts=std::vector<double>();
	if(fractionLifted==0){
		for (int i = 0; i < numberOfNeighborsLifted(nodeID); ++i) {
			size_t neighborID=getNeighborLiftedVertex(nodeID,i);
			//liftedCosts[neighborID]=0;
			liftedCosts.push_back(0);
			liftedIDs.push_back(neighborID);
			liftedIDToOrder[neighborID]=i;
		}
	}
	else{
		for (int i = 0; i < numberOfNeighborsLifted(nodeID); ++i) {
			size_t edgeID=getNeighborLiftedEdge(nodeID,i);
			size_t neighborID=getNeighborLiftedVertex(nodeID,i);
			double cost=ldpInstance.getLiftedEdgeScore(edgeID);
			liftedCosts.push_back(fractionLifted*cost);
			liftedIDs.push_back(neighborID);
			liftedIDToOrder[neighborID]=i;
		}
	}
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::initBaseCosts(double fractionBase){
	baseCosts=std::vector<double>();
	if(fractionBase==0){
		for (int i = 0; i < numberOfNeighborsBase(nodeID); ++i) {
			size_t neighborID=getNeighborBaseVertex(nodeID,i);
			baseCosts.push_back(0);
			baseIDs.push_back(neighborID);
			baseIDToOrder[neighborID]=i;
			//solutionCosts.push_back(0);
//			baseCosts[neighborID]=0;
//			solutionCosts[neighborID]=0;
		}
	}
	else{
		for (int i = 0; i < numberOfNeighborsBase(nodeID); ++i) {
			size_t edgeID=getNeighborBaseEdge(nodeID,i);
			size_t neighborID=getNeighborBaseVertex(nodeID,i);
			double cost=ldpInstance.getEdgeScore(edgeID);
			//baseCosts[neighborID]=fractionBase*cost;
			baseCosts.push_back(fractionBase*cost);
			baseIDs.push_back(neighborID);
			baseIDToOrder[neighborID]=i;

		}
		solutionCosts=std::vector<double>(baseCosts.size()+1);
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
	if(debug())std::cout<<"update values run"<<std::endl;

	//StrForUpdateValues strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID);
	strForUpdateValues.vertexIDStructure.clear();

	strForUpdateValues.valuesStructure.clear();
	strForUpdateValues.setUseAllVertices(true);

	updateValues(strForUpdateValues);

	//std::cout<<"values updated"<<std::endl;

	optimalSolutionBase=strForUpdateValues.optBase;
	optValue=solutionCosts[optimalSolutionBase];
	//std::cout<<"opt base: "<<optimalSolutionBase<<std::endl;

	optimalSolutionLifted=getOptLiftedFromIndexStr(strForUpdateValues);
	//std::cout<<std::endl;

	optLiftedUpToDate=true;
	optBaseUpToDate=true;
}

//TODO make this separate from standard update
template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateValues(StrForUpdateValues& myStr) const{
	updateValues(myStr,nodeNotActiveForStructures);
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateValues(StrForUpdateValues& myStr,size_t vertexIDToIgnore) const{

	//std::cout<<"update values in node "<<nodeID<<std::endl;
	std::unordered_set<size_t> closedVertices;

	//TODO fill valuesStructure with base and lifted endpoints!


	bool lastLayerSet=false;
	size_t lastLayer=0;
	if(vertexIDToIgnore!=nodeNotActiveForStructures){
		lastLayerSet=true;
		lastLayer=ldpInstance.getGroupIndex(vertexIDToIgnore);
	}


	if(lastLayerSet){
		for(auto it=myStr.valuesStructure.begin();it!=myStr.valuesStructure.end();it++){
			if(isInGivenRange(it->first,lastLayer)) it->second=0;
		}
	}
	for (int i = 0; i < liftedIDs.size(); ++i) {
		if(!lastLayerSet||isInGivenRange(liftedIDs[i],lastLayer)){
			myStr.valuesStructure[liftedIDs[i]]=liftedCosts[i];
		}
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
			size_t minValueVertexID=getVertexToReach();

			//std::cout<<"descendants: ";
			for (int i = 0; i < numberOfNeighborsBase(currentNode); ++i) {
				size_t desc=getNeighborBaseVertex(currentNode,i);
				if(desc==vertexIDToIgnore||desc==getVertexToReach()) continue;
				//std::cout<<desc;
				if(isInThisFactorRange(desc)&&myStr.useVertex(desc)){
					if(closedVertices.count(desc)>0||(lastLayerSet&&!isInGivenRange(desc,lastLayer))){  //descendant closed
						//std::cout<<"-cl";
						if(descClosed){
							auto it=myStr.valuesStructure.find(desc);				//	std::cout<<nodeID<<"->"<<bestVertex<<std::endl;find(desc);
							if(it!=myStr.valuesStructure.end()){
								double value=it->second;
								if(minValue>value){
									minValue=value;
									minValueVertexID=desc;
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
					size_t bestIndex=nodeNotActive;

					myStr.solutionCosts[nodeNotActive]=0;
					for (int i = 0; i < baseCosts.size(); ++i) {


						double baseCost=baseCosts[i];
						double valueToAdd=0;
						auto vsIt=myStr.valuesStructure.find(baseIDs[i]);
						if(vsIt!=myStr.valuesStructure.end()){
							valueToAdd=vsIt->second;
						}

						double value=baseCost+nodeCost+valueToAdd;

						myStr.solutionCosts[i]=value;
						if(value<bestValue){
							bestValue=value;
							bestIndex=i;

							//std::cout<<"best vertex"<<std::endl;
						}
						//std::cout<<"end for"<<std::endl;liftedCosts
					}

					myStr.optBase=bestIndex;

					myStr.optValue=bestValue;
				}
				else{

					double valueToStore=minValue;

					myStr.vertexIDStructure[currentNode]=minValueVertexID;
					if(debug())std::cout<<currentNode<<"->"<<minValueVertexID<<", value: "<<valueToStore<<std::endl;

					if((myStr.valuesStructure.count(currentNode)>0)||valueToStore<0){  //store only negative values or values needed to correct solutionCosts
						myStr.valuesStructure[currentNode]+=minValue;

					}

					closedVertices.insert(currentNode); //marking the node as closed.
				}
				nodeStack.pop();

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





template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::oneLiftedMinMarginal(size_t indexOfLiftedEdge)const{
	updateOptimal();

	assert(indexOfLiftedEdge<liftedCosts.size());
	bool isOptimal=false;
	for(size_t optVertex:optimalSolutionLifted){
		if(optVertex==liftedIDs[indexOfLiftedEdge]){
			isOptimal=true;
			break;
		}
	}



	//std::cout<<"lilfted marginal "<<nodeID<<" "<<vertexOfLiftedEdge<<std::endl;

	if(isOptimal){
		//double lbBefore=LowerBound();
		//std::cout<<"OPTIMAL"<<std::endl;


		std::vector<double> localSolutionCosts;
		StrForUpdateValues myStr(baseCosts,liftedCosts,localSolutionCosts,nodeID);
		myStr.copyFromOther(strForUpdateValues);

		myStr.setUseAllVertices(false);

		updateValues(myStr,liftedIDs[indexOfLiftedEdge]);
		//size_t restrictedOptimalSolution=myStr.indexStructure[nodeID];
		double restrictedOptValue=myStr.optValue;

		double valueToReturn=optValue-restrictedOptValue;
		if(valueToReturn>1e-8){
			std::cout<<"wrong min marginal"<<nodeID<<", "<<indexOfLiftedEdge<<", opt vertex but positive value "<<valueToReturn<<std::endl;
		}
		if(debug()) std::cout<<"optimal "<<valueToReturn<<std::endl;


		//std::cout<<"marginal "<<valueToReturn<<"node cost "<<nodeCost<<std::endl;
		return valueToReturn;



	}
	else{


//		for(auto pair:strForUpdateValues.liftedCosts){
//			size_t v=pair.first;
//			double value=liftedCosts.at(v);
//			assert(std::abs(value-pair.second)<eps);
//		}

		std::unordered_map<size_t,size_t> indexStr;
		std::unordered_map<size_t,double> message=bottomUpUpdate(strForUpdateValues,liftedIDs[indexOfLiftedEdge],indexStr);

		auto it =message.begin();
		double messValue=it->second;

		if(debug()) std::cout<<"not optimal "<<messValue<<std::endl;

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
//		for (int i = 0; i < liftedIDs.size(); ++i) {
//			size_t liftedVertexID=liftedIDs[i];
//			double liftedVertexCost=liftedCosts[i];
//			(*pBUValuesStr)[liftedVertexID]=liftedVertexCost;
//		}
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
				auto baseIt=baseIDToOrder.find(currentVertex);
				if(baseIt!=baseIDToOrder.end()){
					double bCost=baseCosts[baseIt->second];
					bestValue=bCost+nodeCost;
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
				auto liftedIt=liftedIDToOrder.find(currentVertex);
				//Assume that open vertex has an entry iff it is lifted
				indexStr[currentVertex]=bestIndex;
				if(liftedIt!=liftedIDToOrder.end()){
					if(onlyOne&&currentVertex!=vertex){
						bestValue+=liftedCosts[liftedIt->second];
					}
					else{
						bestValue+=liftedIt->second;
						double topDownValue=0;
						auto bestTdIt=myStr.vertexIDStructure.find(currentVertex);
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
//
//template<class LDP_INSTANCE>
//inline std::unordered_map<size_t,std::unordered_set<size_t>> ldp_single_node_cut_factor<LDP_INSTANCE>::createCandidateGraph(const StrForUpdateValues& myStr){
//	std::unordered_map<size_t,std::unordered_set<size_t>> candidateGraph;
//	std::unordered_set<size_t> isClosed;
//	for(auto it=myStr.valuesStructure.begin();it!=myStr.valuesStructure.end();++it){
//		size_t vertex=it->first;
//		if(liftedCosts.count(vertex)==0&&isClosed.count(vertex)==0){
//			std::stack<size_t> nodeStack;
//			nodeStack.push(vertex);
//			while(!nodeStack.empty()){
//				size_t currentVertex=nodeStack.top();
//				bool descClosed=true;
//				std::unordered_set<size_t> onlyLiftedDesc;
//				for (int i = 0; i < numberOfNeighborsBase(currentVertex); ++i) {
//					size_t neighbor=getNeighborBaseVertex(currentVertex,i);
//					if(myStr.valuesStructure.count(neighbor)>0){
//						if(liftedCosts.count(neighbor)==0){
//							if(isClosed.count(neighbor)==0){
//								descClosed=false;
//								nodeStack.push(neighbor);
//
//							}
//							else if(descClosed){
//								std::unordered_set<size_t> &neighborLiftedDesc=candidateGraph[neighbor];
//								onlyLiftedDesc.insert(neighborLiftedDesc.begin(),neighborLiftedDesc.end());
//							}
//						}
//						else if(descClosed){
//							onlyLiftedDesc.insert(neighbor);
//						}
//					}
//				}
//				if(descClosed){
//					candidateGraph[currentVertex]=onlyLiftedDesc;
//					nodeStack.pop();
//					isClosed.insert(currentVertex);
//				}
//			}
//		}
//	}
//
//	std::unordered_map<size_t,std::unordered_set<size_t>> finalCandidateGraph;
//
//	for(auto it=myStr.valuesStructure.begin();it!=myStr.valuesStructure.end();++it){
//		size_t vertex=it->first;
//		if(liftedCosts.count(vertex)>0){
//			for (int i = 0; i < numberOfNeighborsBase(vertex); ++i) {
//				size_t neighbor=getNeighborBaseVertex(vertex,i);
//				if(myStr.valuesStructure.count(neighbor)>0){
//					if(liftedCosts.count(neighbor)>0){
//						finalCandidateGraph[vertex].insert(neighbor);
//					}
//					else{
//						std::unordered_set<size_t>& neighborsDesc=candidateGraph[neighbor];
//						finalCandidateGraph[vertex].insert(neighborsDesc.begin(),neighborsDesc.end());
//					}
//				}
//			}
//		}
//	}
//
////		for(auto it=myStr.valuesStructure.begin();it!=myStr.valuesStructure.end();++it){
////			size_t vertex=it->first;
////			//strForUpdateValues.relevantVertices.insert(vertex);
////
////			for (int i = 0; i < numberOfNeighborsBase(vertex); ++i) {
////				size_t neighbor=getNeighborBaseVertex(vertex,i);
////				if(myStr.valuesStructure.count(neighbor)>0){
////					candidateGraph[vertex].insert(neighbor);
////				}
////			}
////		}
//
//		return finalCandidateGraph;
//}

template<class LDP_INSTANCE>
inline std::vector<double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllLiftedMinMarginals(){


	updateOptimal();
	std::unordered_map<size_t,double> liftedMessages;

	double currentOptValue=strForUpdateValues.optValue;

	std::list<size_t> isNotZeroInOpt=optimalSolutionLifted;
	std::unordered_set<size_t> isOneInOpt(isNotZeroInOpt.begin(),isNotZeroInOpt.end());

	std::vector<double> localSolutionCosts=solutionCosts;
	std::vector<double> localLiftedCosts=liftedCosts;
	std::vector<double> localBaseCosts=baseCosts;

	StrForUpdateValues myStr(localBaseCosts,localLiftedCosts,localSolutionCosts,nodeID);
	myStr.copyFromOther(strForUpdateValues);

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
	for(size_t vertexID:liftedIDs){
		if(closedVertices.count(vertexID)==0){

			std::unordered_map<size_t,double> newMessages=bottomUpUpdate(myStr,vertexID,indexStr,&closedVertices,&buValuesStructure);
			liftedMessages.insert(newMessages.begin(),newMessages.end());
		}
	}


	std::vector<double> messagesToOutput=std::vector<double>(liftedCosts.size());
	for (int i = 0; i < messagesToOutput.size(); ++i) {
		messagesToOutput[i]=liftedMessages[liftedIDs[i]];
	}
	return messagesToOutput;
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
