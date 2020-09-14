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
template<class T>
struct ShiftedVector{
public:
    ShiftedVector<T>(size_t boundary1,size_t boundary2,const T& value): //inclusive both min and max vertex
    minVertex(std::min(boundary1,boundary2)),maxVertex(std::max(boundary1,boundary2))
    {
        //const size_t min = std::min(boundary1,boundary2);
        //const size_t max = std::max(boundary1,boundary2);
        //const size_t size = max - min + 1;
        //myVector = new T[size];
        myVector=std::vector<T>(maxVertex-minVertex+1,value);
    }
    ShiftedVector<T>(size_t boundary1,size_t boundary2): //inclusive both min and max vertex
    minVertex(std::min(boundary1,boundary2)),maxVertex(std::max(boundary1,boundary2))
    {
        myVector=std::vector<T>(maxVertex-minVertex+1);
    }
    ShiftedVector<T>(){
        minVertex=0;
        maxVertex=0;
        ShiftedVector(minVertex,maxVertex);
    }


    T& operator[](size_t idx){
        if(debug()&&( idx<minVertex||idx>maxVertex)) std::cout<<"out of bounds, index "<<idx<<", interval "<<minVertex<<","<<maxVertex<<std::endl;
        assert(idx>=minVertex&&idx<=maxVertex);
        size_t shiftedIndex=idx-minVertex;
        return myVector[shiftedIndex];
    }

    void setValue(size_t index,T value){
        assert(index>=minVertex&&index<=maxVertex);
        myVector[index-minVertex]=value;

    }

    T getValue(size_t idx){
        assert(idx>=minVertex&&idx<=maxVertex);
        size_t shiftedIndex=idx-minVertex;
        return myVector[shiftedIndex];
    }

    const T& operator [](size_t idx) const {
        if(debug()&&( idx<minVertex||idx>maxVertex)) std::cout<<"out of bounds, index "<<idx<<", interval "<<minVertex<<","<<maxVertex<<std::endl;
        assert(idx>=minVertex&&idx<=maxVertex);
        return myVector[idx-minVertex];
    }

    void fillWith(const T& value){
        myVector=std::vector<T>(maxVertex-minVertex+1,value);
    }

    void fillWith(const T& value,size_t index1,size_t index2){ //inclusive both boundary indices
        size_t minIndex=std::min(index1,index2);
        size_t maxIndex=std::max(index1,index2);
        assert(minIndex>=minVertex&&maxIndex<=maxVertex);
        for (size_t i = minIndex; i <= maxIndex; ++i) {
            myVector[i-minVertex]=value;
        }
    }

private:
    //T* myVector; // TODO: might be faster depending on whether initializing myVector is expensive
    std::vector<T> myVector;
    size_t minVertex;
    size_t maxVertex;
};



struct StrForUpdateValues{
    std::vector<double> solutionCosts;
    ShiftedVector<double> valuesStructure;
    ShiftedVector<size_t> vertexIDStructure;
//	std::unordered_map<size_t,double> valuesStructure;
//	std::unordered_map<size_t,size_t> vertexIDStructure;
	size_t optBase;
	const std::vector<double>& baseCosts;
	const std::vector<double>& liftedCosts;



	//std::unordered_set<size_t> relevantVertices;
//	bool useAllVertices;
	double optValue;
	const size_t nodeID;
	//size_t optimalSolution;
    //StrForUpdateValues(const std::vector<double>& bCosts,const std::vector<double>& lCosts,std::vector<double>& solCosts,const size_t centralNode):
    //StrForUpdateValues(const std::vector<double>& bCosts,const std::vector<double>& lCosts,const std::vector<double>& solCosts,const size_t centralNode,const size_t boundaryValue,const size_t vertexToReach):
    StrForUpdateValues(const std::vector<double>& bCosts,const std::vector<double>& lCosts,const size_t centralNode,const size_t boundaryValue,const size_t vertexToReach):
	baseCosts(bCosts),
    liftedCosts(lCosts),
    solutionCosts(bCosts.size()+1),
    //solutionCosts(solCosts),
	optValue(0),
    nodeID(centralNode),
    valuesStructure(centralNode,boundaryValue,0),
    vertexIDStructure(centralNode,boundaryValue,vertexToReach),
     optBase(0) //maybe set to node not active
	{

	}


//	bool useVertex(size_t vertex){
//		if(useAllVertices) return true;
//	    return valuesStructure.count(vertex)>0;
//	}

//	void setUseAllVertices(bool value){
//		useAllVertices=value;
//	}
	double getOptimalValue(){
		return optValue;
	}

//	void copyFromOther(StrForUpdateValues& str){
//		vertexIDStructure=str.vertexIDStructure;
//		valuesStructure=str.valuesStructure;
//        //solutionCosts=str.solutionCosts;
//        optValue=str.optValue;
//	}


};

template<class LDP_INSTANCE>
class ldp_single_node_cut_factor
{
public:
	//constexpr static std::size_t nodeNotActive = std::numeric_limits<std::size_t>::max();

	//template<class LPD_STRUCT> ldp_single_node_cut_factor(const LPD_STRUCT& ldpStruct);
	ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut);
	ldp_single_node_cut_factor(const ldp_single_node_cut_factor&);
    ~ldp_single_node_cut_factor(){
       // std::cout<<"destructor "<<nodeID<<std::endl;
    }

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
        assert(primalBase_<baseIDs.size());
        return baseIDs.at(primalBase_);
	}

	const std::unordered_set<size_t>& getPrimalLiftedIndices() const {
		return primalLifted_;
	}

	const bool isActiveInPrimalLifted(size_t vertex) const {
		return primalLifted_.count(vertex)>0;
		//return primalLifted_.at(vertex);
	}

    const size_t getLiftedIDToOrder(size_t vertexID) const{
        auto it=liftedIDToOrder.find(vertexID);
        assert(it!=liftedIDToOrder.end());
        return it->second;
	}

    const size_t getBaseIDToOrder(size_t vertexID) const{
        auto it=baseIDToOrder.find(vertexID);
        assert(it!=baseIDToOrder.end());
        return it->second;
    }
    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(baseCosts, liftedCosts); }

    //auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }//TODO change this. This will not work with so many variables
    auto export_variables() { return std::tie(baseCosts, liftedCosts); }

	void updateCostSimple(const double value,const size_t vertexIndex,bool isLifted);
	void updateNodeCost(const double value);
    double getOneBaseEdgeMinMarginal(size_t vertex) const;
    std::vector<double> getAllBaseMinMarginals()const;

    std::vector<double> getAllLiftedMinMarginals(std::vector<double>* pLocalBaseCosts=nullptr) const;

	double getNodeMinMarginal()const;
	double oneLiftedMinMarginal(size_t indexOfLiftedEdge)const;

	const std::vector<size_t>& getBaseIDs() const {
		return baseIDs;
	}

    size_t getBaseID(size_t order) const {
        return baseIDs.at(order);
    }

	const std::vector<size_t>& getLiftedIDs() const {
		return liftedIDs;
	}

    size_t getLiftedID(size_t order) const {
        return liftedIDs.at(order);
    }

	size_t getNodeNotActive() const {
		return nodeNotActive;
	}

	const std::size_t nodeID;
	//const size_t nodeNotActive;
	size_t primalBase_;
	std::unordered_set<size_t> primalLifted_;

	//TODO use this
	//std::vector<bool> primalLifted_;


    void print() const{
        std::cout<<nodeID<<":";
        if(isOutFlow){
            std::cout<<"out";
        }
        else{
            std::cout<<"in";
        }
    }



    const LDP_INSTANCE& getLdpInstance() const{
        return ldpInstance;
    }



private:
  //  void updateValues() const;
    void updateValues(StrForUpdateValues& myStr, const size_t vertexIDToIgnore) const;
	void updateValues(StrForUpdateValues& myStr) const;
    std::unordered_map<size_t,double> bottomUpUpdate(const StrForUpdateValues& myStr, const size_t vertex, ShiftedVector<size_t>& indexStr, ShiftedVector<bool>* pClosedVert=nullptr, ShiftedVector<double>* pBUValuesStr=nullptr)const;
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
//	bool isInThisFactorRange(const size_t vertexID) const {
//		assert(vertexID < baseGraph.numberOfVertices());
//		if(isOutFlow){
//			//if(debug()) std::cout<<"vertex: "<<vertexID<<", max l. "<<maxLayer<<", layer of vertex "<<ldpInstance.getGroupIndex(vertexID);
//			if(vertexID==ldpInstance.getTerminalNode()) return true;
//			else return ldpInstance.getGroupIndex(vertexID)<=maxLayer;
//		}
//		else{
//			if(vertexID==ldpInstance.getSourceNode()) return true;
//			else return ldpInstance.getGroupIndex(vertexID)>=minLayer;
//		}
//	}

//	bool isInGivenRange(const size_t nodeIndex,const size_t boundLayer) const {
//		assert(nodeIndex < baseGraph.numberOfVertices());
//		if(isOutFlow){
//			return ldpInstance.getGroupIndex(nodeIndex)<=boundLayer;
//		}
//		else{
//			return ldpInstance.getGroupIndex(nodeIndex)>=boundLayer;
//		}
//    }


    bool isInGivenInterval(const size_t nodeIndex,const size_t boundaryIndex) const {
        assert(nodeIndex < baseGraph.numberOfVertices());
        if(isOutFlow){
            return nodeIndex<=boundaryIndex;
        }
        else{
            return nodeIndex>=boundaryIndex;
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
	size_t getNeighborLiftedEdge(size_t firstNode,size_t neighborIndex)const {
		if(isOutFlow){
			return liftedGraph.edgeFromVertex(firstNode,neighborIndex);
		}
		else{
			return liftedGraph.edgeToVertex(firstNode,neighborIndex);
		}
	}

    size_t getNeighborLiftedVertex(size_t firstNode,size_t neighborIndex)const {
		if(isOutFlow){
			return liftedGraph.vertexFromVertex(firstNode,neighborIndex);
		}
		else{
			return liftedGraph.vertexToVertex(firstNode,neighborIndex);
		}
	}

	size_t numberOfNeighborsLifted(size_t nodeIndex)const {
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


	std::pair<bool,size_t> findEdgeBase(size_t firstNode,size_t secondNode)const {
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

    size_t getMoreDistantNode(const size_t& v1,const size_t& v2)const {
        if(isOutFlow){
            return std::max(v1,v2);
        }
        else{
            return std::min(v1,v2);
        }
    }


    mutable size_t optimalSolutionBase; //Contains vertices w.r.t. their order for central vertex
    mutable std::list<size_t> optimalSolutionLifted;  //Contains IDs of vertices, not their order w.r.t. central vertex!

	//std::set<size_t> optimalSolutionLifted;

	//std::set<size_t,decltype(vertexCompare)> mySet;

	std::size_t minLayer;
	std::size_t maxLayer;

    std::size_t mostDistantNeighborID;

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

    mutable bool solutionCostsUpToDate;
    mutable bool optValueUpToDate;

	mutable double optValue;

//	mutable StrForUpdateValues strForUpdateValues;


	std::vector<size_t> baseIDs;
	std::vector<size_t> liftedIDs;

	size_t nodeNotActive;


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
mostDistantNeighborID(sncFactor.mostDistantNeighborID),
//strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID),
nodeNotActive(sncFactor.nodeNotActive),
solutionCostsUpToDate(sncFactor.solutionCostsUpToDate),
 optValueUpToDate(sncFactor.optValueUpToDate)
{
  //  std::cout<<"copy constructor "<<nodeID<<std::endl;

	primalBase_=sncFactor.primalBase_;
	primalLifted_=sncFactor.primalLifted_;

	optimalSolutionBase=sncFactor.optimalSolutionBase;
	optimalSolutionLifted=sncFactor.optimalSolutionLifted;

	minLayer=sncFactor.minLayer;
	maxLayer=sncFactor.maxLayer;

	nodeCost=sncFactor.nodeCost;
	optValue=sncFactor.optValue;

//	strForUpdateValues.copyFromOther(sncFactor.strForUpdateValues);


	optLiftedUpToDate=sncFactor.optLiftedUpToDate;
	optBaseUpToDate=sncFactor.optBaseUpToDate;

    nodeNotActiveForStructures=sncFactor.nodeNotActiveForStructures;
	baseIDs=sncFactor.baseIDs;
	liftedIDs=sncFactor.liftedIDs;
	liftedIDToOrder=sncFactor.liftedIDToOrder;
	baseIDToOrder=sncFactor.baseIDToOrder;


}


template<class LDP_INSTANCE>
inline  ldp_single_node_cut_factor<LDP_INSTANCE>::ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut):
baseGraph(ldpInst.getGraph()),
liftedGraph(ldpInst.getGraphLifted()),
nodeID(nID),
ldpInstance(ldpInst),
isOutFlow(isOut)
//strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID)
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
    mostDistantNeighborID=nodeID;
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
    optValueUpToDate=false;
    solutionCostsUpToDate=false;

//	if(debug()){
//			std::cout<<"factor "<<nodeID<<" created."<<std::endl;
//		}



}


template<class LDP_INSTANCE>
inline std::list<size_t> ldp_single_node_cut_factor<LDP_INSTANCE>::getOptLiftedFromIndexStr(const StrForUpdateValues& myStr) const{

	std::list<size_t> optLifted;
    //if(debug())	std::cout<<"opt lifted: "<<nodeID<<std::endl;
	double optValueComputed=0;
	if(myStr.optBase!=nodeNotActive){
		optValueComputed=myStr.baseCosts.at(myStr.optBase);
		optValueComputed+=nodeCost;



        size_t vertexInOptimalPath=baseIDs.at(myStr.optBase);
        //size_t vertexInOptimalPath=myStr.vertexIDStructure.at(baseIDs[myStr.optBase]);
        bool hasOptDescendant=vertexInOptimalPath!=getVertexToReach();
		while(hasOptDescendant){
        //    std::cout<<"opt vertex "<<vertexInOptimalPath<<std::endl;

			if(liftedIDToOrder.count(vertexInOptimalPath)>0){

				optLifted.push_back(vertexInOptimalPath);
				if(debug()){

                    double toAdd=myStr.liftedCosts.at(liftedIDToOrder.at(vertexInOptimalPath));
                  //  std::cout<<"is lifted, value "<<toAdd<<std::endl;
					optValueComputed+=toAdd;
				}
			}
            size_t vert=myStr.vertexIDStructure[vertexInOptimalPath];
            hasOptDescendant=vert!=getVertexToReach();
			if(hasOptDescendant){
                vertexInOptimalPath=vert;
			}
		}
	}


	return optLifted;

}




//template<class LDP_INSTANCE>
//inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateOptimal() const{
//	if(!optLiftedUpToDate){

//		//if(debug())std::cout<<"updating lifted str."<<std::endl;

//		updateValues();
//	}
//	else if(!optBaseUpToDate){
//		//if(debug())std::cout<<"updating base str."<<std::endl;
//		optimalSolutionBase=nodeNotActive;
//		double minValue=0;
//		for (int i = 0; i < solutionCosts.size()-1; ++i) {
//			double value=solutionCosts[i];
//			if(value<minValue){
//				minValue=value;
//                optimalSolutionBase=i;

//			}
//		}


//        StrForUpdateValues strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID,mostDistantNeighborID);

//		strForUpdateValues.optBase=optimalSolutionBase;
//		strForUpdateValues.optValue=minValue;
//		optimalSolutionLifted=getOptLiftedFromIndexStr(strForUpdateValues);
//		optBaseUpToDate=true;
//        optValue=strForUpdateValues.optValue;
//        solutionCosts=strForUpdateValues.solutionCosts;
//	}
//    //optValue=strForUpdateValues.optValue;
//}



template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setPrimalLifted(std::unordered_set<size_t>& verticesOfActiveEdges) {
	primalLifted_=verticesOfActiveEdges;
}





template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::EvaluatePrimal() const{
    //double value=0;
    //if(debug())
    //std::cout<<"evaluate primal "<<nodeID<<": ";

//    if(nodeID==16&&isOutFlow){
//        std::cout<<"base costs "<<nodeID<<std::endl;
//        for (int i = 0; i < baseCosts.size(); ++i) {
//            std::cout<<i<<": "<<baseCosts[i]<<std::endl;
//        }
//        std::cout<<"lifted costs "<<std::endl;
//        for (int j = 0; j < liftedCosts.size(); ++j) {
//            std::cout<<j<<": "<<liftedCosts[j]<<std::endl;
//        }
//        std::cout<<"primal solution "<<primalBase_<<std::endl;
//        std::cout<<" lifted solution< ";
//    }
    if(primalBase_==nodeNotActive){
      //  std::cout<<0.0<<std::endl;
        return 0;
    }
    else{
        double value=nodeCost;
        value+=baseCosts.at(primalBase_);

        for(size_t node:primalLifted_){
            value+=liftedCosts.at(node);
        //    if(nodeID==16&&isOutFlow) std::cout<<node<<", ";
        }
        //if(nodeID==16&&isOutFlow) std::cout<<std::endl;
      //  std::cout<<value<<std::endl;
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
inline void ldp_single_node_cut_factor<LDP_INSTANCE>:: updateOptimal()const {
    if(!optValueUpToDate){
        if(solutionCostsUpToDate){
            optValue=solutionCosts[0];
            for(size_t i=0;i<solutionCosts.size();i++){
                if(optValue>solutionCosts[i]){
                    optValue=solutionCosts[i];
                    optimalSolutionBase=i;
                }
            }
        }
        else{
            StrForUpdateValues myStr(baseCosts,liftedCosts,nodeID,mostDistantNeighborID,getVertexToReach());
            updateValues(myStr);
            optValue=myStr.optValue;
            solutionCosts=myStr.solutionCosts;
            optimalSolutionBase=myStr.optBase;
          //  std::cout<<"updated structures "<<std::endl;

        }
    }
    optValueUpToDate=true;
    solutionCostsUpToDate=true;
  //  std::cout<<"finished update optimal"<<std::endl;

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
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getOneBaseEdgeMinMarginal(const size_t index)const{
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


template<class LDP_INSTANCE>
inline std::vector<double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllBaseMinMarginals() const{
	updateOptimal();
	//if(debug()) std::cout<<"output min marginals"<<std::endl;
	std::vector<double> minMarginals(baseCosts.size());
	if(optimalSolutionBase==nodeNotActive){
		for (int i = 0; i < solutionCosts.size()-1; ++i) {
			minMarginals[i]=solutionCosts[i];
			//if(debug()) std::cout<<i<<" "<<minMarginals[i]<<std::endl;

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
			//if(debug()) std::cout<<i<<" "<<minMarginals[i]<<std::endl;
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
        optBaseUpToDate=false; //TODO sometimes can stay true
        optValueUpToDate=false;
	}
	else{ //update in lifted edge

		assert(vertexIndex<liftedCosts.size());
		liftedCosts[vertexIndex]+=value;

		optLiftedUpToDate=false;
		optBaseUpToDate=false;
        optValueUpToDate=false;
        solutionCostsUpToDate=false;


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
	liftedIDs=std::vector<size_t>();
	liftedIDToOrder.clear();
	if(fractionLifted==0){
		for (int i = 0; i < numberOfNeighborsLifted(nodeID); ++i) {
			size_t neighborID=getNeighborLiftedVertex(nodeID,i);
            if(neighborID!=getVertexToReach()){
                mostDistantNeighborID=getMoreDistantNode(mostDistantNeighborID,neighborID);
            }
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
            if(neighborID!=getVertexToReach()){
                mostDistantNeighborID=getMoreDistantNode(mostDistantNeighborID,neighborID);
            }
			double cost=ldpInstance.getLiftedEdgeScore(edgeID);
			liftedCosts.push_back(fractionLifted*cost);
			liftedIDs.push_back(neighborID);
			liftedIDToOrder[neighborID]=i;
		}
	}
    optLiftedUpToDate=false;
    optValueUpToDate=false;
    solutionCostsUpToDate=false;
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::initBaseCosts(double fractionBase){
	baseCosts=std::vector<double>();
	baseIDs=std::vector<size_t>();
	baseIDToOrder.clear();
	if(fractionBase==0){
		for (int i = 0; i < numberOfNeighborsBase(nodeID); ++i) {
			size_t neighborID=getNeighborBaseVertex(nodeID,i);
            if(neighborID!=getVertexToReach()){
                mostDistantNeighborID=getMoreDistantNode(mostDistantNeighborID,neighborID);
            }
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
            if(neighborID!=getVertexToReach()){
                mostDistantNeighborID=getMoreDistantNode(mostDistantNeighborID,neighborID);
            }
			double cost=ldpInstance.getEdgeScore(edgeID);
			//baseCosts[neighborID]=fractionBase*cost;
			baseCosts.push_back(fractionBase*cost);
			baseIDs.push_back(neighborID);
			baseIDToOrder[neighborID]=i;

		}

	}
    solutionCosts=std::vector<double>(baseCosts.size()+1);
    //	for (int i = 0; i < numberOfNeighborsLifted(nodeID); ++i) {
	//		size_t edgeID=getNeighborLiftedEdge(nodeID,i);
	//		size_t neighborID=getNeighborLiftedVertex(nodeID,i);
	//		double cost=ldpInstance.getLiftedEdgeScore(edgeID);
	//		liftedCosts[neighborID]=fractionLifted*cost;
	//	}
    optBaseUpToDate=false;
    optValueUpToDate=false;
    solutionCostsUpToDate=false;

}


template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::LowerBound() const{//TODO store info about how valuesStructures changed. At least max time layer of changed lifted edge
    //if(debug()) std::cout<<"lower bound "<<nodeID<<std::endl;
    updateOptimal();
    return optValue;

}



//template<class LDP_INSTANCE>
//inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateValues() const{
//	//if(debug())std::cout<<"update values run"<<nodeID<<std::endl;

//    StrForUpdateValues strForUpdateValues(baseCosts,liftedCosts,solutionCosts,nodeID,mostDistantNeighborID);
//    //strForUpdateValues.vertexIDStructure.clear();

//    //strForUpdateValues.valuesStructure.clear();
//	strForUpdateValues.setUseAllVertices(true);

//	updateValues(strForUpdateValues);

//	//std::cout<<"values updated"<<std::endl;

//	optimalSolutionBase=strForUpdateValues.optBase;
//    optValue=solutionCosts.at(optimalSolutionBase);
//	//std::cout<<"opt base: "<<optimalSolutionBase<<std::endl;

//	optimalSolutionLifted=getOptLiftedFromIndexStr(strForUpdateValues);
//	//std::cout<<std::endl;

//	optLiftedUpToDate=true;
//	optBaseUpToDate=true;
//}

//TODO make this separate from standard update
template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateValues(StrForUpdateValues& myStr) const{
	updateValues(myStr,nodeNotActiveForStructures);
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateValues(StrForUpdateValues& myStr,const size_t vertexIDToIgnore) const{

    //std::cout<<"update values in node "<<nodeID<<std::endl;
   // std::cout<<"most distant "<<mostDistantNeighborID<<std::endl;
    //std::cout<<"to reach "<<getVertexToReach()<<std::endl;
    ShiftedVector<bool> closedVertices(nodeID,mostDistantNeighborID);


    bool vertexToIgnoreSet=false;
    size_t lastVertex=mostDistantNeighborID;
	if(vertexIDToIgnore!=nodeNotActiveForStructures){
        vertexToIgnoreSet=true;
        lastVertex=vertexIDToIgnore;
	}


    if(vertexToIgnoreSet){
        myStr.valuesStructure.fillWith(0,nodeID,vertexIDToIgnore);
    }
//	else{  //This should not be needed. Now done in str. constructor
//        myStr.valuesStructure.fillWith(0);
//        myStr.vertexIDStructure.fillWith(getVertexToReach());
//	//	std::cout<<"clear vs "<<std::endl;
//	}
    for (int i = 0; i < liftedIDs.size(); ++i) {
        //assert(liftedIDs.at(i) < baseGraph.numberOfVertices());
        if(!vertexToIgnoreSet||isInGivenInterval(liftedIDs.at(i),lastVertex)){

            myStr.valuesStructure[liftedIDs.at(i)]=myStr.liftedCosts.at(i);

		}
	}


	std::stack<size_t> nodeStack;
	nodeStack.push(nodeID);

    //if(debug()) std::cout<<"update values before while "<<std::endl;
	while(!nodeStack.empty()){
		size_t currentNode=nodeStack.top();
        // if(debug()) std::cout<<"current node "<<currentNode<<std::endl;
        if(closedVertices.getValue(currentNode)){
			nodeStack.pop();
		}
		else{

			bool descClosed=true;
			double minValue=0;
			//std::unordered_set<size_t> minValueIndices;
			size_t minValueVertexID=getVertexToReach();

			//std::cout<<"descendants: ";
			for (int i = 0; i < numberOfNeighborsBase(currentNode); ++i) {
				size_t desc=getNeighborBaseVertex(currentNode,i);
                if(desc==vertexIDToIgnore||desc==getVertexToReach()) continue;

                if(isInGivenInterval(desc,mostDistantNeighborID)){

                    assert(desc<baseGraph.numberOfVertices());
                    if(closedVertices.getValue(desc)||(vertexToIgnoreSet&&!isInGivenInterval(desc, vertexIDToIgnore))){  //descendant closed
                        if(descClosed){
                            double value=myStr.valuesStructure[desc];
                            if(minValue>value){
                                minValue=value;
                                minValueVertexID=desc;
                            }
                        }
                    }
					else {  //descendant not closed

						nodeStack.push(desc);
                       // std::cout<<"push "<<desc<<std::endl;

						descClosed=false;

					}
                }


            }
			if(descClosed){ //Close node if all descendants are closed


               // std::cout<<"desc closed"<<std::endl;
				if(currentNode==nodeID){  //all nodes closed, compute solution values
					double bestValue=0;
					size_t bestIndex=nodeNotActive;

					myStr.solutionCosts[nodeNotActive]=0;
                    for (size_t i = 0; i < myStr.baseCosts.size(); ++i) {

                       // std::cout<<"base cost for "<<i<<", "<<baseIDs[i]<<", vertex to reach "<<getVertexToReach()<<std::endl;

                        double baseCost=myStr.baseCosts[i];

                        double valueToAdd=0;
                        if(baseIDs[i]!=getVertexToReach()){
                            valueToAdd=myStr.valuesStructure[baseIDs[i]];
                        }
						double value=baseCost+nodeCost+valueToAdd;

                        myStr.solutionCosts.at(i)=value;
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

                    assert(currentNode<baseGraph.numberOfVertices());
					myStr.vertexIDStructure[currentNode]=minValueVertexID;
                    myStr.valuesStructure[currentNode]+=minValue;

                    closedVertices.setValue(currentNode,true); //marking the node as closed.
                    //std::cout<<"closed value "<<closedVertices.getValue(currentNode)<<std::endl;
				}
				nodeStack.pop();

			}
		}

	}
 //   std::cout<<"update values finished"<<std::endl;



}



template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::oneLiftedMinMarginal(size_t indexOfLiftedEdge)const{
    assert(indexOfLiftedEdge<liftedCosts.size());


    StrForUpdateValues strForUpdateValues(baseCosts,liftedCosts,nodeID,mostDistantNeighborID,getVertexToReach());
    updateValues(strForUpdateValues);
    double origOptValue=strForUpdateValues.optValue;

    std::list<size_t> optimalSolutionLifted=getOptLiftedFromIndexStr(strForUpdateValues);

	bool isOptimal=false;
	for(size_t optVertex:optimalSolutionLifted){
        if(optVertex==liftedIDs.at(indexOfLiftedEdge)){
			isOptimal=true;
			//if(debug()) std::cout<<"is optimal"<<std::endl;
			break;
		}
	}



  //  if(debug())std::cout<<"lifted marginal "<<nodeID<<" "<<indexOfLiftedEdge<<", optimal "<<isOptimal<<std::endl;

	if(isOptimal){

        updateValues(strForUpdateValues,liftedIDs.at(indexOfLiftedEdge));

        //size_t restrictedOptimalSolution=myStr.indexStructure[nodeID];
        double restrictedOptValue=strForUpdateValues.optValue;

        double valueToReturn=origOptValue-restrictedOptValue;
		if(debug()&&valueToReturn>1e-8){
			std::cout<<"wrong min marginal"<<nodeID<<", "<<indexOfLiftedEdge<<", opt vertex but positive value "<<valueToReturn<<std::endl;
		}


		return valueToReturn;



	}
	else{

        ShiftedVector<size_t> indexStr(nodeID,mostDistantNeighborID,getVertexToReach());

        std::unordered_map<size_t,double> message=bottomUpUpdate(strForUpdateValues,liftedIDs[indexOfLiftedEdge],indexStr);

		auto it =message.begin();
		double messValue=it->second;

		if(debug()){
			if(messValue<-eps){
				std::cout<<"wrong min marginal "<<nodeID<<", index "<<indexOfLiftedEdge<<". Not optimal, value "<<messValue<<std::endl;
			}
		}

        return messValue;
	}

}


template<class LDP_INSTANCE>
inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::bottomUpUpdate(const StrForUpdateValues& myStr,const size_t vertex,ShiftedVector<size_t>& indexStr,ShiftedVector<bool>* pClosedVert,ShiftedVector<double>* pBUValuesStr)const{
  //  std::cout<<"bottom up update from "<<vertex<<std::endl;
    bool onlyOne=pClosedVert==nullptr;


	std::unordered_map<size_t,double> messages;
	if(onlyOne){
        pClosedVert=new ShiftedVector<bool> (nodeID,mostDistantNeighborID);
        pBUValuesStr=new ShiftedVector<double>(nodeID,mostDistantNeighborID);
	}

    ShiftedVector<bool>& closedVertices=*pClosedVert;
    ShiftedVector<double>& buValuesStr=*pBUValuesStr;



	std::stack<size_t> myStack;
	myStack.push(vertex);

	while(!myStack.empty()){
		size_t currentVertex=myStack.top();
        if(closedVertices.getValue(currentVertex)){
			myStack.pop();
		}
		else{
			bool predClosed=true;
			for (int i = 0; i < numberOfNeighborsBaseRev(currentVertex); ++i) {
				size_t pred=getNeighborBaseVertexRev(currentVertex,i);
				if(pred==nodeID) continue;
                if(reachable(nodeID,pred)&&!closedVertices.getValue(pred)){
					predClosed=false;
					myStack.push(pred);
				}
			}
			if(predClosed){
				size_t bestIndex=nodeID;
				double bestValue=std::numeric_limits<double>::infinity();
				auto baseIt=baseIDToOrder.find(currentVertex);
				if(baseIt!=baseIDToOrder.end()){
                    double bCost=myStr.baseCosts.at(baseIt->second);
					bestValue=bCost+nodeCost;
				}
				for (int i = 0; i < numberOfNeighborsBaseRev(currentVertex); ++i) {
                    size_t pred=getNeighborBaseVertexRev(currentVertex,i);
                    if(pred==nodeID||!reachable(nodeID,pred)) continue;
                    assert(closedVertices.getValue(pred));
                    double value=buValuesStr[pred];
                    if(value<bestValue){
                        bestValue=value;
                        bestIndex=pred;  //TODO check
                    }

                }
                auto liftedIt=liftedIDToOrder.find(currentVertex);
				//Assume that open vertex has an entry iff it is lifted

				indexStr[currentVertex]=bestIndex;
				if(liftedIt!=liftedIDToOrder.end()){
					if(onlyOne&&currentVertex!=vertex){
                        bestValue+=myStr.liftedCosts.at(liftedIt->second);
					}
					else{
//                        bestValue+=myStr.liftedCosts.at(liftedIt->second);
//						double topDownValue=0;

//                        size_t bestTd=myStr.vertexIDStructure[currentVertex];
//						if(bestTd!=getVertexToReach()){
//                            topDownValue=myStr.valuesStructure[bestTd];
//						}

//						double restrictedOpt=topDownValue+bestValue;
//						double delta=restrictedOpt-myStr.optValue;

//						messages[currentVertex]=delta;
//                      //  std::cout<<"message "<<currentVertex<<": "<<delta<<std::endl;
//						bestValue-=delta;


                        double topDownValueOfDesc=0;

                        size_t bestDesc=myStr.vertexIDStructure[currentVertex];
                        if(bestDesc!=getVertexToReach()){
                            topDownValueOfDesc=myStr.valuesStructure[bestDesc];
                        }
                        double topDownCurrent=myStr.valuesStructure[currentVertex];

                        double restrictedOpt=topDownCurrent+bestValue;
                        double delta=restrictedOpt-myStr.optValue;
                        bestValue=myStr.optValue-topDownValueOfDesc;

                        messages[currentVertex]=delta;
                      //  std::cout<<"message "<<currentVertex<<": "<<delta<<std::endl;


					}
				}

                closedVertices.setValue(currentVertex,true);
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
inline std::vector<double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllLiftedMinMarginals(std::vector<double>* pLocalBaseCosts) const{



   // std::cout<<"is out "<<isOutFlow<<std::endl;
	std::unordered_map<size_t,double> liftedMessages;


//	std::vector<double> localSolutionCosts=solutionCosts;
	std::vector<double> localLiftedCosts=liftedCosts;
    std::vector<double> localBaseCosts;
    if(pLocalBaseCosts==nullptr){
        localBaseCosts=baseCosts;
    }
    else{
        localBaseCosts=*pLocalBaseCosts;
    }

    StrForUpdateValues myStr(localBaseCosts,localLiftedCosts,nodeID,mostDistantNeighborID,getVertexToReach());
    updateValues(myStr);
    //std::cout<<"further all lifted min marginals "<<std::endl;

    std::list<size_t> isNotZeroInOpt=getOptLiftedFromIndexStr(myStr);
    std::unordered_set<size_t> isOneInOpt(isNotZeroInOpt.begin(),isNotZeroInOpt.end());


    double origOptValue=myStr.optValue;
    double minMarginalsImproving=0;
    double currentOptValue=myStr.optValue;

   // std::cout<<"orig opt value "<<origOptValue<<std::endl;
    //myStr.setUseAllVertices(false);


    auto listIt=isNotZeroInOpt.begin();
	while(!isNotZeroInOpt.empty()){
		size_t vertexToClose=*listIt;
//        if(debug()) std::cout<<"process vertex "<<vertexToClose<<std::endl;

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
        size_t orderToClose=liftedIDToOrder.at(vertexToClose);
        localLiftedCosts[orderToClose]-=delta;
		liftedMessages[vertexToClose]=delta;
        minMarginalsImproving+=delta;
		currentOptValue=newOpt;
      //  std::cout<<"curren opt "<<currentOptValue<<std::endl;

	}

//	myStr.setUseAllVertices(true);

   // std::cout<<"min marginals improve "<<minMarginalsImproving<<std::endl;
    myStr.valuesStructure.fillWith(0);  //TODO implement this in reset in StrForUpdate
    myStr.vertexIDStructure.fillWith(getVertexToReach());
	updateValues(myStr);

    assert(currentOptValue+minMarginalsImproving-origOptValue>-eps);

   // std::cout<<"computed opt "<<myStr.optValue<<std::endl;
    assert(std::abs(currentOptValue-myStr.optValue)<eps);
    ShiftedVector<double> buValuesStructure(nodeID,mostDistantNeighborID,std::numeric_limits<double>::max());
    ShiftedVector<bool> closedVertices(nodeID,mostDistantNeighborID,false);

    for(size_t optVertex:isOneInOpt){
        assert(optVertex<baseGraph.numberOfVertices());
        buValuesStructure[optVertex]=currentOptValue-myStr.valuesStructure[optVertex]+localLiftedCosts[liftedIDToOrder.at(optVertex)];
        closedVertices.setValue(optVertex,true);
	}
    ShiftedVector<size_t> indexStr(nodeID,mostDistantNeighborID,nodeID);
	for(size_t vertexID:liftedIDs){
        if(!closedVertices.getValue(vertexID)){
			std::unordered_map<size_t,double> newMessages=bottomUpUpdate(myStr,vertexID,indexStr,&closedVertices,&buValuesStructure);
			liftedMessages.insert(newMessages.begin(),newMessages.end());
            if(debug()){
                for(auto pair:newMessages){
                    size_t v=pair.first;
                    double val=pair.second;
                    localLiftedCosts.at(liftedIDToOrder.at(v))-=val;
                }
            }
		}
	}


	std::vector<double> messagesToOutput=std::vector<double>(liftedCosts.size());
	for (int i = 0; i < messagesToOutput.size(); ++i) {
		messagesToOutput[i]=liftedMessages[liftedIDs[i]];
	}

    if(debug()){
        StrForUpdateValues myStr2(localBaseCosts,localLiftedCosts,nodeID,mostDistantNeighborID,getVertexToReach());

        updateValues(myStr2);
        assert(std::abs(myStr2.optValue-currentOptValue)<eps);
        assert(myStr2.optValue+minMarginalsImproving-origOptValue>-eps);
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
        if(debug()) std::cout<<"repam left "<<l.nodeID<<": "<<l.getLiftedID(right_node)<<":"<<msg<<std::endl;
		assert(msg_dim == 0);
		l.updateCostSimple(msg,right_node,true);
	}

	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
	{
        if(debug()) std::cout<<"repam right "<<r.nodeID<<": "<<r.getLiftedID(left_node)<<":"<<msg<<std::endl;
		assert(msg_dim == 0);
		r.updateCostSimple(msg,left_node,true);
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
	{
        if(debug()) std::cout<<"message to left "<<r.nodeID<<": "<<r.getLiftedID(left_node)<<std::endl;
		const double delta = r.oneLiftedMinMarginal(left_node);
		msg[0] -= omega * delta;
    //    if(debug()) std::cout<<"sent "<<r.nodeID<<": "<<r.getLiftedID(left_node)<<std::endl;
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_right(const SINGLE_NODE_CUT_FACTOR& l, MSG& msg, const double omega)
	{
        if(debug()) std::cout<<"message to right "<<l.nodeID<<": "<<l.getLiftedID(right_node)<<std::endl;
		const double delta = l.oneLiftedMinMarginal(right_node);
		msg[0] -= omega * delta;
	}



    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG_ARRAY>
    static void SendMessagesToLeft(const SINGLE_NODE_CUT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
    {
        if(debug()){
            std::cout<<"running get all lifted marginals to left "<<r.nodeID<<std::endl;
        }
        const std::vector<double> msg_vec = r.getAllLiftedMinMarginals();
//        if(debug()){
//            std::cout<<"obtained all lifted marginals, size "<<msg_vec.size()<<std::endl;
//        }
        for(auto it=msg_begin; it!=msg_end; ++it)
        {
            auto& msg = (*it).GetMessageOp();
            const size_t left_node = msg.left_node;
            const size_t right_node = msg.right_node;
            (*it)[0] -= omega * msg_vec.at(left_node);
        }
//        if(debug()){
//            std::cout<<"messages added "<<std::endl;
//        }
    }

    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG_ARRAY>
    static void SendMessagesToRight(const SINGLE_NODE_CUT_FACTOR& l, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
    {
        if(debug()){
            std::cout<<"running get all lifted marginals to right "<<l.nodeID<<std::endl;
        }
        const std::vector<double> msg_vec = l.getAllLiftedMinMarginals();
//        if(debug()){
//            std::cout<<"obtained all lifted marginals, size "<<msg_vec.size()<<std::endl;
//        }
        for(auto it=msg_begin; it!=msg_end; ++it)
        {
            auto& msg = (*it).GetMessageOp();
            const size_t left_node = msg.left_node;
            const size_t right_node = msg.right_node;
//            if(debug()){
//               // auto& origVal=(*it)[0];
//                std::cout<<"source node"<<l.nodeID<<"to node "<<l.getLiftedID(right_node)<<", to subtract "<<msg_vec.at(right_node)<<std::endl;
//            }
            //(*it)[0] -= omega * msg_vec.at(right_node);
            (*it).operator[](0)-= omega * msg_vec.at(right_node);
        }
        if(debug()){
            std::cout<<"messages added "<<std::endl;
        }
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
