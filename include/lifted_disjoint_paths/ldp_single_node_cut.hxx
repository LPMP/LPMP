#pragma once
#include <andres/graph/digraph.hxx>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <array>
#include <list>
#include <set>
#include <config.hxx>
#include "ldp_directed_graph.hxx"
#include<lifted_disjoint_paths/ldp_functions.hxx>

namespace LPMP {



template<class INSTANCE>
struct StrForTopDownUpdate{

    //Structure storing results of topDownUpdate
    StrForTopDownUpdate(const std::vector<double>& bCosts,const std::vector<double>& lCosts,const INSTANCE&instance, const size_t centralNodeID,const size_t mostDistantNodeID,const size_t vertexToReach):
	baseCosts(bCosts),
    liftedCosts(lCosts),
    solutionCosts(bCosts.size()+1),
	optValue(0),
    nodeID(centralNodeID)

	{
        size_t first=std::min(centralNodeID,mostDistantNodeID);
        size_t last=std::max(centralNodeID,mostDistantNodeID)+1;


        fillWithValue<size_t>(instance.sncNeighborStructure,first,last,vertexToReach);
        fillWithValue<double>(instance.sncTDStructure,first,last,0);


        optBaseIndex=bCosts.size();
	}

    const size_t nodeID;
    const std::vector<double>& baseCosts;
    const std::vector<double>& liftedCosts;

    std::vector<double> solutionCosts;
    size_t optBaseIndex;
    double optValue;


};

template<class LDP_INSTANCE>
class ldp_single_node_cut_factor
{
public:
    ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut);

    //Getting and setting lifted and base costs
    void initBaseCosts(double fractionBase);

	void initLiftedCosts(double fractionLifted);

	const std::vector<double>& getLiftedCosts() const {
		return liftedCosts;
    }

	const std::vector<double>& getBaseCosts() const {
		return baseCosts;
	}

    //Methods for getting and setting primal solution
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


    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(baseCosts, liftedCosts); }

    auto export_variables() { return std::tie(baseCosts, liftedCosts); }

    //Updating costs
    void updateEdgeCost(const double value,const size_t vertexIndex,bool isLifted);
	void updateNodeCost(const double value);


    //Computing lower bound and min marginals
    double LowerBound() const;

    double getOneBaseEdgeMinMarginal(size_t vertex) const;
    double getOneLiftedMinMarginal(size_t indexOfLiftedEdge)const;
    double getNodeMinMarginal()const;

    std::vector<double> getAllBaseMinMarginals()const;
    std::vector<double> getAllBaseMinMarginalsForMCF() const;
    std::vector<double> getAllLiftedMinMarginals(const std::vector<double> *pLocalBaseCosts=nullptr) const;


    //Accessing maps between local and global node indices
    const std::vector<size_t>& getBaseIDs() const {
		return baseIDs;
	}

    size_t getBaseID(size_t index) const {
        return baseIDs.at(index);
    }

	const std::vector<size_t>& getLiftedIDs() const {
		return liftedIDs;
	}

    size_t getLiftedID(size_t index) const {
        return liftedIDs.at(index);
    }

    const size_t getLiftedIDToOrder(size_t vertexID) const{
        auto it=liftedIDToOrder.find(vertexID);
        assert(it!=liftedIDToOrder.end());
        return it->second;
    }

    const size_t getBaseIDToOrder(size_t vertexID) const{
        auto it=baseIDToIndex.find(vertexID);
        assert(it!=baseIDToIndex.end());
        return it->second;
    }

    //Printing node information for debugging purposes
    void print() const{
        std::cout<<nodeID<<":";
        if(isOutFlow){
            std::cout<<"out";
        }
        else{
            std::cout<<"in";
        }
    }

         const std::size_t nodeID;


private:
    //Methods used in computing lower bound or min marginals
    void topDownUpdate(StrForTopDownUpdate<LDP_INSTANCE>& myStr, const size_t vertexIDToIgnore) const;
    void topDownUpdate(StrForTopDownUpdate<LDP_INSTANCE>& myStr) const;
    std::unordered_map<size_t,double> bottomUpUpdate(const StrForTopDownUpdate<LDP_INSTANCE>& myStr, const size_t vertex, const ShiftedVector<char> &verticesInScope, bool onlyOne)const;
    void updateOptimal() const;

    //Obtain IDs of vertices of lifted edges that are part of a found optimal solution
    std::list<size_t> getOptLiftedFromIndexStr(const StrForTopDownUpdate<LDP_INSTANCE>& myStr)const;


    //Methods for exploring graph structures
//    size_t getNeighborBaseVertex(size_t firstNode,size_t neighborIndex) const{
//        assert(firstNode < baseGraph.numberOfVertices());
//        if(isOutFlow){
//            return baseGraph.vertexFromVertex(firstNode,neighborIndex);
//        }
//        else{
//            return baseGraph.vertexToVertex(firstNode,neighborIndex);
//        }
//    }
//    size_t numberOfNeighborsBase(const size_t nodeIndex) const {
//        assert(nodeIndex < baseGraph.numberOfVertices());
//        if(isOutFlow){
//            return baseGraph.numberOfEdgesFromVertex(nodeIndex);
//        }
//        else{
//            return baseGraph.numberOfEdgesToVertex(nodeIndex);
//        }
//    }
//    size_t numberOfNeighborsBaseRev(const size_t nodeIndex) const {
//        assert(nodeIndex < baseGraph.numberOfVertices());
//        if(!isOutFlow){
//            return baseGraph.numberOfEdgesFromVertex(nodeIndex);
//        }
//        else{
//            return baseGraph.numberOfEdgesToVertex(nodeIndex);
//        }
//    }


    bool isInGivenInterval(const size_t nodeIndex,const size_t boundaryIndex) const {
        assert(nodeIndex < ldpInstance.getNumberOfVertices());
        if(isOutFlow){
            return nodeIndex<=boundaryIndex;
        }
        else{
            return nodeIndex>=boundaryIndex;
        }
    }

//    size_t getNeighborBaseEdge(size_t firstNode,size_t neighborIndex)const{
//        if(isOutFlow){
//            return baseGraph.edgeFromVertex(firstNode,neighborIndex);
//        }
//        else{
//            return baseGraph.edgeToVertex(firstNode,neighborIndex);
//        }
//    }

//    size_t getNeighborBaseVertexRev(size_t firstNode,size_t neighborIndex)const{
//        if(!isOutFlow){
//            return baseGraph.vertexFromVertex(firstNode,neighborIndex);
//        }
//        else{
//            return baseGraph.vertexToVertex(firstNode,neighborIndex);
//        }
//    }

    const size_t* neighborsBegin(const size_t& nodeIndex)const{
        if(isOutFlow){
            return ldpBaseGraph.forwardNeighborsBegin(nodeIndex);
        }
        else {
            return ldpBaseGraph.backwardNeighborsBegin(nodeIndex);
        }
    }

    const size_t* neighborsEnd(const size_t& nodeIndex)const{
        if(isOutFlow){
            return ldpBaseGraph.forwardNeighborsEnd(nodeIndex);
        }
        else {
            return ldpBaseGraph.backwardNeighborsEnd(nodeIndex);
        }
    }

    const size_t* neighborsRevBegin(const size_t& nodeIndex)const{
        if(!isOutFlow){
            return ldpBaseGraph.forwardNeighborsBegin(nodeIndex);
        }
        else {
            return ldpBaseGraph.backwardNeighborsBegin(nodeIndex);
        }
    }

    const size_t* neighborsRevEnd(const size_t& nodeIndex)const{
        if(!isOutFlow){
            return ldpBaseGraph.forwardNeighborsEnd(nodeIndex);
        }
        else {
            return ldpBaseGraph.backwardNeighborsEnd(nodeIndex);
        }
    }

//	size_t getNeighborLiftedEdge(size_t firstNode,size_t neighborIndex)const {
//		if(isOutFlow){
//			return liftedGraph.edgeFromVertex(firstNode,neighborIndex);
//		}
//		else{
//			return liftedGraph.edgeToVertex(firstNode,neighborIndex);
//		}
//	}

//    size_t getNeighborLiftedVertex(size_t firstNode,size_t neighborIndex)const {
//		if(isOutFlow){
//			return liftedGraph.vertexFromVertex(firstNode,neighborIndex);
//		}
//		else{
//			return liftedGraph.vertexToVertex(firstNode,neighborIndex);
//		}
//	}

//	size_t numberOfNeighborsLifted(size_t nodeIndex)const {
//		if(isOutFlow){
//			return liftedGraph.numberOfEdgesFromVertex(nodeIndex);
//		}
//		else{
//			return liftedGraph.numberOfEdgesToVertex(nodeIndex);
//		}
//	}

	bool reachable(size_t firstVertex,size_t secondVertex)const{
		if(isOutFlow){
			return ldpInstance.isReachable(firstVertex,secondVertex);
		}
		else{
			return ldpInstance.isReachable(secondVertex,firstVertex);
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

     bool isLiftedVertex(size_t v)const{
         if(isOutFlow){
             return ldpInstance.existLiftedEdge(nodeID,v);
         }
         else{
             return ldpInstance.existLiftedEdge(v,nodeID);
         }
     }
     const bool isOutFlow; //is it outgoing flow
     size_t minVertex;
     size_t maxVertex;

     //References to graph structures
    // const andres::graph::Digraph<>& baseGraph;
    // const andres::graph::Digraph<>& liftedGraph;
     const LDP_INSTANCE& ldpInstance;
     const LdpDirectedGraph& ldpBaseGraph;
     const LdpDirectedGraph& ldpLiftedGraph;

     //costs
     std::vector<double> baseCosts;
     std::vector<double> liftedCosts;
     double nodeCost;

     //Variables storing primal solutions
     size_t primalBase_;
     std::unordered_set<size_t> primalLifted_;

     //ID of the node most distant from the central node that is contained in this factor
     std::size_t mostDistantNeighborID;
     //Index in solutionCosts for storing node not active value
     size_t nodeNotActive;


     //Best values achievable by activating given base edges, last value corresponds to no edge active
     mutable std::vector<double> solutionCosts;
     //Index of the current optimal base edge or equal to nodeNotActive
     mutable size_t optBaseIndex;
     //Current optimal value
     mutable double optValue;

     //Is vector solutionCosts up to date
     mutable bool solutionCostsUpToDate;
     //Is optValue up to date
     mutable bool optValueUpToDate;

     //Maps from indices in this SNC factor to global node IDs and back
     std::vector<size_t> baseIDs;
     std::vector<size_t> liftedIDs;
     std::unordered_map<size_t,size_t> baseIDToIndex;
     std::unordered_map<size_t,size_t> liftedIDToOrder;


};


template<class LDP_INSTANCE>
inline  ldp_single_node_cut_factor<LDP_INSTANCE>::ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut):
//baseGraph(ldpInst.getGraph()),
//liftedGraph(ldpInst.getGraphLifted()),
nodeID(nID),
ldpInstance(ldpInst),
isOutFlow(isOut),
ldpBaseGraph(ldpInst.getMyGraph()),
ldpLiftedGraph(ldpInst.getMyGraphLifted())

{


    mostDistantNeighborID=nodeID;
    initBaseCosts(0);
    initLiftedCosts(0);

    if(isOutFlow){
        minVertex=nodeID;
        maxVertex=mostDistantNeighborID;
    }
    else{
        maxVertex=nodeID;
        minVertex=mostDistantNeighborID;
    }

	nodeNotActive=baseCosts.size();

	primalBase_=nodeNotActive;  //corresponds to no edge active
    optBaseIndex=nodeNotActive;
	nodeCost=0;
	optValue=0;

    optValueUpToDate=false;
    solutionCostsUpToDate=false;



}


template<class LDP_INSTANCE>
inline std::list<size_t> ldp_single_node_cut_factor<LDP_INSTANCE>::getOptLiftedFromIndexStr(const StrForTopDownUpdate<LDP_INSTANCE>& myStr) const{

	std::list<size_t> optLifted;
    double optValueComputed=0;
    if(myStr.optBaseIndex!=nodeNotActive){
        optValueComputed=myStr.baseCosts.at(myStr.optBaseIndex);
        optValueComputed+=nodeCost;

        size_t vertexInOptimalPath=baseIDs.at(myStr.optBaseIndex);

        while(vertexInOptimalPath!=getVertexToReach()){

            if(isLiftedVertex(vertexInOptimalPath)){

                optLifted.push_back(vertexInOptimalPath);
                double toAdd=myStr.liftedCosts.at(liftedIDToOrder.at(vertexInOptimalPath));
                optValueComputed+=toAdd;

            }
           // vertexInOptimalPath=myStr.topDownVertexIDStructure[vertexInOptimalPath];
             vertexInOptimalPath=ldpInstance.sncNeighborStructure[vertexInOptimalPath];
		}
        assert(std::abs(optValueComputed-myStr.optValue)<eps);
    }

    return optLifted;

}





template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setPrimalLifted(std::unordered_set<size_t>& verticesOfActiveEdges) {
	primalLifted_=verticesOfActiveEdges;
}





template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::EvaluatePrimal() const{
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
inline void ldp_single_node_cut_factor<LDP_INSTANCE>:: updateOptimal()const {
    if(!optValueUpToDate){
        if(solutionCostsUpToDate){
            optValue=solutionCosts[0];
            optBaseIndex=0;
            for(size_t i=0;i<solutionCosts.size();i++){
                if(optValue>solutionCosts[i]){
                    optValue=solutionCosts[i];
                    optBaseIndex=i;
                }
            }
        }
        else{
            StrForTopDownUpdate myStr(baseCosts,liftedCosts,ldpInstance,nodeID,mostDistantNeighborID,getVertexToReach());
            topDownUpdate(myStr);
            optValue=myStr.optValue;
            solutionCosts=myStr.solutionCosts;
            optBaseIndex=myStr.optBaseIndex;
            solutionCostsUpToDate=true;
          //  std::cout<<"updated structures "<<std::endl;

        }
        optValueUpToDate=true;
    }


  //  std::cout<<"finished update optimal"<<std::endl;

}

template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getNodeMinMarginal()const{
	//std::cout<<"node min marginal "<<nodeID<<" "<<isOutFlow<<std::endl;
	updateOptimal();
    if(optBaseIndex!=nodeNotActive){
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
    if(optBaseIndex!=index){
        return solutionCosts[index]-solutionCosts[optBaseIndex];
	}
	else{
		double secondBest=std::numeric_limits<double>::max();
		for (int i = 0; i < solutionCosts.size(); ++i) {
            if(i==optBaseIndex) continue;
			if(solutionCosts[i]<secondBest){
				secondBest=solutionCosts[i];
			}
		}
		return solutionCosts[index]-secondBest;
	}
}



template<class LDP_INSTANCE>
inline std::vector<double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllBaseMinMarginalsForMCF() const{
    updateOptimal();
    auto it=solutionCosts.end();
    it--;
    std::vector<double> toReturn(solutionCosts.begin(),it);
    return toReturn;

}



template<class LDP_INSTANCE>
inline std::vector<double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllBaseMinMarginals() const{
	updateOptimal();
	//if(debug()) std::cout<<"output min marginals"<<std::endl;
	std::vector<double> minMarginals(baseCosts.size());
    if(optBaseIndex==nodeNotActive){
		for (int i = 0; i < solutionCosts.size()-1; ++i) {
			minMarginals[i]=solutionCosts[i];
			//if(debug()) std::cout<<i<<" "<<minMarginals[i]<<std::endl;

		}
	}
	else{
		double secondBest=std::numeric_limits<double>::infinity();
        double optValue=solutionCosts.at(optBaseIndex);
		for (int i = 0; i < solutionCosts.size(); ++i) {
            if(i==optBaseIndex) continue;
			if(solutionCosts[i]<secondBest){
				secondBest=solutionCosts[i];
			}

		}
		for (int i = 0; i < solutionCosts.size()-1; ++i) {
			minMarginals[i]=solutionCosts[i]-secondBest;

		}
    }
	return minMarginals;
}

template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateNodeCost(const double value){
	nodeCost+=value;
	for (int i = 0; i < solutionCosts.size()-1; ++i) {
		(solutionCosts[i]+=value);
	}
    optValueUpToDate=false;

}

template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateEdgeCost(const double value,const size_t vertexIndex,bool isLifted){//Only cost change
	if(!isLifted){ //update in base edge
		assert(vertexIndex<baseCosts.size());
		baseCosts[vertexIndex]+=value;
		solutionCosts[vertexIndex]+=value;
        optValueUpToDate=false;
	}
	else{ //update in lifted edge

		assert(vertexIndex<liftedCosts.size());
		liftedCosts[vertexIndex]+=value;

        optValueUpToDate=false;
        solutionCostsUpToDate=false;


	}
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::initLiftedCosts(double fractionLifted){
	liftedCosts=std::vector<double>();
	liftedIDs=std::vector<size_t>();
	liftedIDToOrder.clear();
    const LdpDirectedGraph& myLiftedGraph=ldpInstance.getMyGraphLifted();

    if(isOutFlow){
        const double* costIt=myLiftedGraph.forwardCostBegin(nodeID);
        const size_t* edgeIt=myLiftedGraph.forwardNeighborsBegin(nodeID);
        size_t counter=0;
        for (;edgeIt!=myLiftedGraph.forwardNeighborsEnd(nodeID);costIt++, edgeIt++) {
            assert(costIt!=myLiftedGraph.forwardCostEnd(nodeID));
            liftedCosts.push_back((*costIt)*fractionLifted);
            liftedIDs.push_back(*edgeIt);
            liftedIDToOrder[*edgeIt]=counter;
            if(*edgeIt!=ldpInstance.getTerminalNode()){
                mostDistantNeighborID=std::max(mostDistantNeighborID,*edgeIt);
            }
            counter++;
        }
    }
    else{
        const double* costIt=myLiftedGraph.backwardCostBegin(nodeID);
        const size_t* edgeIt=myLiftedGraph.backwardNeighborsBegin(nodeID);
        size_t counter=0;
        for (;edgeIt!=myLiftedGraph.backwardNeighborsEnd(nodeID);costIt++, edgeIt++) {
            assert(costIt!=myLiftedGraph.backwardCostEnd(nodeID));
            liftedCosts.push_back((*costIt)*fractionLifted);
            liftedIDs.push_back(*edgeIt);
            liftedIDToOrder[*edgeIt]=counter;
            if(*edgeIt!=ldpInstance.getSourceNode()){
                mostDistantNeighborID=std::min(mostDistantNeighborID,*edgeIt);
            }
            counter++;
        }
    }

//	if(fractionLifted==0){
//		for (int i = 0; i < numberOfNeighborsLifted(nodeID); ++i) {
//			size_t neighborID=getNeighborLiftedVertex(nodeID,i);
//            if(neighborID!=getVertexToReach()){
//                mostDistantNeighborID=getMoreDistantNode(mostDistantNeighborID,neighborID);
//            }
//			//liftedCosts[neighborID]=0;
//			liftedCosts.push_back(0);
//			liftedIDs.push_back(neighborID);
//			liftedIDToOrder[neighborID]=i;
//		}
//	}
//	else{
//		for (int i = 0; i < numberOfNeighborsLifted(nodeID); ++i) {
//			size_t edgeID=getNeighborLiftedEdge(nodeID,i);
//			size_t neighborID=getNeighborLiftedVertex(nodeID,i);
//            if(neighborID!=getVertexToReach()){
//                mostDistantNeighborID=getMoreDistantNode(mostDistantNeighborID,neighborID);
//            }
//			double cost=ldpInstance.getLiftedEdgeScore(edgeID);
//			liftedCosts.push_back(fractionLifted*cost);
//			liftedIDs.push_back(neighborID);
//			liftedIDToOrder[neighborID]=i;
//		}
//	}

    optValueUpToDate=false;
    solutionCostsUpToDate=false;
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::initBaseCosts(double fractionBase){
	baseCosts=std::vector<double>();
	baseIDs=std::vector<size_t>();
    baseIDToIndex.clear();


    const LdpDirectedGraph& myBaseGraph=ldpInstance.getMyGraph();

    if(isOutFlow){
        const double* costIt=myBaseGraph.forwardCostBegin(nodeID);
        const size_t* edgeIt=myBaseGraph.forwardNeighborsBegin(nodeID);
        size_t counter=0;
        for (;edgeIt!=myBaseGraph.forwardNeighborsEnd(nodeID);costIt++, edgeIt++) {
            assert(costIt!=myBaseGraph.forwardCostEnd(nodeID));
            baseCosts.push_back((*costIt)*fractionBase);
            baseIDs.push_back(*edgeIt);
            baseIDToIndex[*edgeIt]=counter;
            if(*edgeIt!=ldpInstance.getTerminalNode()){
                mostDistantNeighborID=std::max(mostDistantNeighborID,*edgeIt);
            }
            counter++;
        }
    }
    else{
        const double* costIt=myBaseGraph.backwardCostBegin(nodeID);
        const size_t* edgeIt=myBaseGraph.backwardNeighborsBegin(nodeID);
        size_t counter=0;
        for (;edgeIt!=myBaseGraph.backwardNeighborsEnd(nodeID);costIt++, edgeIt++) {
            assert(costIt!=myBaseGraph.backwardCostEnd(nodeID));
            baseCosts.push_back((*costIt)*fractionBase);
            baseIDs.push_back(*edgeIt);
            baseIDToIndex[*edgeIt]=counter;
            if(*edgeIt!=ldpInstance.getSourceNode()){
                mostDistantNeighborID=std::min(mostDistantNeighborID,*edgeIt);
            }
            counter++;
        }
    }




//	if(fractionBase==0){
//		for (int i = 0; i < numberOfNeighborsBase(nodeID); ++i) {
//			size_t neighborID=getNeighborBaseVertex(nodeID,i);
//            if(neighborID!=getVertexToReach()){
//                mostDistantNeighborID=getMoreDistantNode(mostDistantNeighborID,neighborID);
//            }
//			baseCosts.push_back(0);
//			baseIDs.push_back(neighborID);
//            baseIDToIndex[neighborID]=i;
//		}
//	}
//	else{
//		for (int i = 0; i < numberOfNeighborsBase(nodeID); ++i) {
//			size_t edgeID=getNeighborBaseEdge(nodeID,i);
//			size_t neighborID=getNeighborBaseVertex(nodeID,i);
//            if(neighborID!=getVertexToReach()){
//                mostDistantNeighborID=getMoreDistantNode(mostDistantNeighborID,neighborID);
//            }
//			double cost=ldpInstance.getEdgeScore(edgeID);
//            baseCosts.push_back(fractionBase*cost);
//			baseIDs.push_back(neighborID);
//            baseIDToIndex[neighborID]=i;

//		}

//	}
    solutionCosts=std::vector<double>(baseCosts.size()+1);

    optValueUpToDate=false;
    solutionCostsUpToDate=false;

}


template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::LowerBound() const{//TODO store info about how valuesStructures changed. At least max time layer of changed lifted edge
    updateOptimal();
    return optValue;

}




//TODO make this separate from standard update
template<class LDP_INSTANCE>
void ldp_single_node_cut_factor<LDP_INSTANCE>::topDownUpdate(StrForTopDownUpdate<LDP_INSTANCE>& myStr) const{
     topDownUpdate(myStr,getVertexToReach());
}


template<class LDP_INSTANCE>
void ldp_single_node_cut_factor<LDP_INSTANCE>::topDownUpdate(StrForTopDownUpdate<LDP_INSTANCE>& myStr,const size_t vertexIDToIgnore) const{

    //ShiftedVector<char> closedVertices(nodeID,mostDistantNeighborID);

    //fillWithValue<char>(ldpInstance.sncClosedVertices,minVertex,maxVertex+1,0);

    bool vertexToIgnoreSet=false;
    size_t lastVertex=mostDistantNeighborID;
    if(vertexIDToIgnore!=getVertexToReach()){
        vertexToIgnoreSet=true;
        lastVertex=vertexIDToIgnore;
        size_t minV=0;
        size_t maxV=0;
        if(isOutFlow){
            minV=vertexIDToIgnore;
            maxV=mostDistantNeighborID+1;
            fillWithValue<char>(ldpInstance.sncClosedVertices,minV,maxV,1);
            fillWithValue<char>(ldpInstance.sncClosedVertices,nodeID,minV,0);
        }
        else{
            minV=mostDistantNeighborID;
            maxV=vertexIDToIgnore+1;
            fillWithValue<char>(ldpInstance.sncClosedVertices,minV,maxV,1);
            fillWithValue<char>(ldpInstance.sncClosedVertices,maxV,nodeID+1,0);

        }

	}
    else{
        fillWithValue<char>(ldpInstance.sncClosedVertices,minVertex,maxVertex+1,0);
    }


    if(vertexToIgnoreSet){
       // myStr.topDownValuesStructure.fillWith(0,nodeID,vertexIDToIgnore);
        fillWithValue<double>(ldpInstance.sncTDStructure,std::min(nodeID,vertexIDToIgnore),std::max(nodeID,vertexIDToIgnore)+1,0);
    }

    //Store all lifted costs to top down values structure
    for (int i = 0; i < liftedIDs.size(); ++i) {
        if(!vertexToIgnoreSet||isInGivenInterval(liftedIDs.at(i),lastVertex)){
            //myStr.topDownValuesStructure[liftedIDs.at(i)]=myStr.liftedCosts.at(i);
            ldpInstance.sncTDStructure[liftedIDs.at(i)]=myStr.liftedCosts.at(i);
		}
	}


	std::stack<size_t> nodeStack;
	nodeStack.push(nodeID);


	while(!nodeStack.empty()){
		size_t currentNode=nodeStack.top();

        if(ldpInstance.sncClosedVertices[currentNode]){
			nodeStack.pop();
		}
		else{

			bool descClosed=true;
            double bestDescValue=0;
            size_t bestDescVertexID=getVertexToReach();

            //Traverse all base neighbors in direction out from central node
            //Not closed vertices store to the stack. Search for best descendant in closed neighbors
            const size_t* vertexIt=neighborsBegin(currentNode);
            const size_t* end=neighborsEnd(currentNode);

         //   const double* costIt=ldpBaseGraph.forwardCostBegin(currentNode);
            for (;vertexIt!=end;vertexIt++) {

               size_t desc=*vertexIt;
//			for (int i = 0; i < numberOfNeighborsBase(currentNode); ++i) {
//				size_t desc=getNeighborBaseVertex(currentNode,i);
                if(desc==vertexIDToIgnore||desc==getVertexToReach()) continue;

                if(isInGivenInterval(desc,mostDistantNeighborID)){

                    //assert(desc<baseGraph.numberOfVertices());
                    if(ldpInstance.sncClosedVertices[desc]){  //descendant closed
                        if(descClosed){
                            //double value=myStr.topDownValuesStructure[desc];
                            double value=ldpInstance.sncTDStructure[desc];
                            if(bestDescValue>value){
                                bestDescValue=value;
                                bestDescVertexID=desc;
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

				if(currentNode==nodeID){  //all nodes closed, compute solution values
                    double bestSolutionValue=0;
                    size_t bestSolutionIndex=nodeNotActive;

					myStr.solutionCosts[nodeNotActive]=0;
                    for (size_t i = 0; i < myStr.baseCosts.size(); ++i) {
                        double baseCost=myStr.baseCosts[i];
                        if(!vertexToIgnoreSet||baseIDs[i]!=vertexIDToIgnore){

                            double valueToAdd=0;
                            if(baseIDs[i]!=getVertexToReach()){
                               // valueToAdd=myStr.topDownValuesStructure[baseIDs[i]];
                                 valueToAdd=ldpInstance.sncTDStructure[baseIDs[i]];
                            }
                            double value=baseCost+nodeCost+valueToAdd;

                            myStr.solutionCosts.at(i)=value;
                            if(value<bestSolutionValue){
                                bestSolutionValue=value;
                                bestSolutionIndex=i;

                            }
                        }

                    }

                    myStr.optBaseIndex=bestSolutionIndex;
                    myStr.optValue=bestSolutionValue;

                }
                else{
                    //Close current node, store best node value and pointer to the best descendant

                    assert(!vertexToIgnoreSet||bestDescVertexID!=vertexIDToIgnore);
                    assert(currentNode<ldpInstance.getNumberOfVertices());

                    //adding the descendant value to the already stored lifted cost of node
//                    myStr.topDownValuesStructure[currentNode]+=bestDescValue;
//                    myStr.topDownVertexIDStructure[currentNode]=bestDescVertexID;
                    ldpInstance.sncTDStructure[currentNode]+=bestDescValue;
                    ldpInstance.sncNeighborStructure[currentNode]=bestDescVertexID;

                    ldpInstance.sncClosedVertices[currentNode]=1; //marking the node as closed.

				}
				nodeStack.pop();

			}
		}

	}
 //   std::cout<<"update values finished"<<std::endl;

   // return closedVertices;

}



template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getOneLiftedMinMarginal(size_t indexOfLiftedEdge)const{
    assert(indexOfLiftedEdge<liftedCosts.size());
    //std::cout<<"One lifted min marginal"<<std::endl;


    StrForTopDownUpdate strForUpdateValues(baseCosts,liftedCosts,ldpInstance,nodeID,mostDistantNeighborID,getVertexToReach());
    topDownUpdate(strForUpdateValues);
    ShiftedVector<char> verticesInScope(minVertex,maxVertex,ldpInstance.sncClosedVertices);
    double origOptValue=strForUpdateValues.optValue;

   // std::cout<<"one lifted min marginal"<<std::endl;
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

        //If is optimal, run topDownUpdate with ignoring this vertex
        topDownUpdate(strForUpdateValues,liftedIDs.at(indexOfLiftedEdge));
        double restrictedOptValue=strForUpdateValues.optValue;

        double valueToReturn=origOptValue-restrictedOptValue;
		if(debug()&&valueToReturn>1e-8){
			std::cout<<"wrong min marginal"<<nodeID<<", "<<indexOfLiftedEdge<<", opt vertex but positive value "<<valueToReturn<<std::endl;
		}
		return valueToReturn;



	}
	else{

        //If it is not optimal, find the best possible solution containing this vertex active via bottomUpUpdate
//        ShiftedVector<size_t> indexStr(nodeID,mostDistantNeighborID,getVertexToReach());
//        ShiftedVector<double> buValuesStructure(nodeID,mostDistantNeighborID,0);
        fillWithValue<size_t>(ldpInstance.sncBUNeighborStructure,minVertex,maxVertex+1,getVertexToReach());
        fillWithValue<double>(ldpInstance.sncBUStructure,minVertex,maxVertex,0);

        for (int i = 0; i < baseIDs.size(); ++i) {
            if(baseIDs.at(i)==getVertexToReach()) continue;
//            indexStr[baseIDs.at(i)]=nodeID;
//            buValuesStructure[baseIDs.at(i)]=baseCosts.at(i);
            ldpInstance.sncBUNeighborStructure[baseIDs.at(i)]=nodeID;
            ldpInstance.sncBUStructure[baseIDs.at(i)]=baseCosts.at(i);
        }

        std::unordered_map<size_t,double> message=bottomUpUpdate(strForUpdateValues,liftedIDs[indexOfLiftedEdge],verticesInScope,true);

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
inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::bottomUpUpdate(const StrForTopDownUpdate<LDP_INSTANCE>& myStr,const size_t vertex,const ShiftedVector<char>& verticesInScope,bool onlyOne)const{
  //  std::cout<<"bottom up update from "<<vertex<<std::endl;
   // bool onlyOne=pClosedVert==nullptr; //If running from getOneLiftedMinMarginal, we want min marginal only for the input vertex

	std::unordered_map<size_t,double> messages;
	if(onlyOne){
        fillWithValue<char>(ldpInstance.sncClosedVertices,minVertex,maxVertex+1,0);
    }

    //ShiftedVector<char>& closedVertices=*pClosedVert;

    //Start DFS from the input vertex in the direction towards the central node
	std::stack<size_t> myStack;
	myStack.push(vertex);

	while(!myStack.empty()){
		size_t currentVertex=myStack.top();
        //if(closedVertices[currentVertex]){
        if(ldpInstance.sncClosedVertices[currentVertex]){
			myStack.pop();
		}
		else{
			bool predClosed=true;
            const size_t* vertexIt=neighborsRevBegin(currentVertex);
            const size_t* end=neighborsRevEnd(currentVertex);
            for (;vertexIt!=end;vertexIt++) {
                 size_t pred=*vertexIt;
//			for (int i = 0; i < numberOfNeighborsBaseRev(currentVertex); ++i) {
//				size_t pred=getNeighborBaseVertexRev(currentVertex,i);
                if(pred==nodeID||!verticesInScope.isWithinBounds(pred)||!verticesInScope[pred]) continue;
                if(!ldpInstance.sncClosedVertices[pred]){
					predClosed=false;
					myStack.push(pred);
				}
			}
            if(predClosed){

                size_t bestIndex=getVertexToReach();
                double bestValue=std::numeric_limits<double>::infinity();
//                if(bottomUpVertexIDStructure[currentVertex]==nodeID){  //This holds for endpoints of base edges
//                    bestValue=bottomUpValuesStructure[currentVertex]+nodeCost;
                if(ldpInstance.sncBUNeighborStructure[currentVertex]==nodeID){  //This holds for endpoints of base edges
                    bestValue=ldpInstance.sncBUStructure[currentVertex]+nodeCost;
                    bestIndex=nodeID;
                }
                const size_t* vertexIt=neighborsRevBegin(currentVertex);
                const size_t* end=neighborsRevEnd(currentVertex);
                for (;vertexIt!=end;vertexIt++) {
                    size_t pred=*vertexIt;

//                for (int i = 0; i < numberOfNeighborsBaseRev(currentVertex); ++i) { //Select the best possible neighbor
//                    size_t pred=getNeighborBaseVertexRev(currentVertex,i);
                    if(pred==nodeID||!verticesInScope.isWithinBounds(pred)||!verticesInScope[pred]) continue;
                    assert(ldpInstance.sncClosedVertices[pred]>0);
                   // double value=bottomUpValuesStructure[pred];
                     double value=ldpInstance.sncBUStructure[pred];
                    if(value<bestValue){
                        bestValue=value;
                        bestIndex=pred;  //TODO check
                    }

                }

                //bottomUpVertexIDStructure[currentVertex]=bestIndex;
                ldpInstance.sncBUNeighborStructure[currentVertex]=bestIndex;
                if(isLiftedVertex(currentVertex)){
                    if(onlyOne&&currentVertex!=vertex){  //For the case of getting one min marginal
                        bestValue+=myStr.liftedCosts.at(liftedIDToOrder.at(currentVertex));
					}
					else{
//
                        double topDownValueOfDesc=0;

                        //size_t bestDesc=myStr.topDownVertexIDStructure[currentVertex]; //Best neighbor in the direction from the central node
                        size_t bestDesc=ldpInstance.sncNeighborStructure[currentVertex]; //Best neighbor in the direction from the central node
                        if(bestDesc!=getVertexToReach()){
                           // topDownValueOfDesc=myStr.topDownValuesStructure[bestDesc];
                             topDownValueOfDesc=ldpInstance.sncTDStructure[bestDesc];
                        }
                        //double topDownCurrent=myStr.topDownValuesStructure[currentVertex];  //Best top down value
                        double topDownCurrent=ldpInstance.sncTDStructure[currentVertex];  //Best top down value

                        double restrictedOpt=topDownCurrent+bestValue;  //Best top down value of currentVertex plus bottom up value of best predecessor
                        double delta=restrictedOpt-myStr.optValue;
                        bestValue=myStr.optValue-topDownValueOfDesc;  //Compute node's bottom up value after changing its lifted cost by delta

                        messages[currentVertex]=delta;
                      //  std::cout<<"message "<<currentVertex<<": "<<delta<<std::endl;


					}
				}

                ldpInstance.sncClosedVertices[currentVertex]=1;
                ldpInstance.sncBUStructure[currentVertex]=bestValue;
				myStack.pop();

			}
		}
	}
//	if(onlyOne){
//		delete(pClosedVert);
//        pClosedVert=nullptr;

//	}
	return messages;


}


template<class LDP_INSTANCE>
inline std::vector<double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllLiftedMinMarginals(const std::vector<double>* pLocalBaseCosts) const{

	std::unordered_map<size_t,double> liftedMessages;


	std::vector<double> localLiftedCosts=liftedCosts;
    std::vector<double> localBaseCosts;
    if(pLocalBaseCosts==nullptr){
        localBaseCosts=baseCosts;
    }
    else{//If other value of base costs needs to be used. E.g. for getting min marginals of both base and lifted edges
        localBaseCosts=*pLocalBaseCosts;
    }

    //First, compute optimal value
    StrForTopDownUpdate myStr(localBaseCosts,localLiftedCosts,ldpInstance,nodeID,mostDistantNeighborID,getVertexToReach());
    topDownUpdate(myStr);
    double origOptValue=myStr.optValue;


    //All vertices that are not zero in any optimal solution
    std::list<size_t> isNotZeroInOpt=getOptLiftedFromIndexStr(myStr);
    //All vertices that are one in at least one of the optimal solutions
    std::unordered_set<size_t> isOneInOpt(isNotZeroInOpt.begin(),isNotZeroInOpt.end());

//    if(myStr.optBaseIndex!=nodeNotActive){
//    std::cout<<"opt base "<< baseIDs.at(myStr.optBaseIndex)<<"opt path "<<std::endl;
//    for(auto it=isNotZeroInOpt.begin();it!=isNotZeroInOpt.end();it++){
//        std::cout<<*it<<",";
//    }
//    std::cout<<std::endl;
//    }
//    else{
//        std::cout<<"opt not active "<<std::endl;
//    }

    double minMarginalsImproving=0;

    double currentOptValue=myStr.optValue;
   // std::cout<<"opt value "<<currentOptValue<<std::endl;
    auto listIt=isNotZeroInOpt.begin();

    //Obtaining min marginals for nodes that are active in all optimal solutions
    while(!isNotZeroInOpt.empty()){
        size_t vertexToClose=*listIt;
       // std::cout<<"closing "<<vertexToClose<<std::endl;

        //Obtaining best solution while ignoring vertexToClose
        topDownUpdate(myStr,vertexToClose);
        double newOpt=myStr.optValue;
        std::list<size_t> secondBest=getOptLiftedFromIndexStr(myStr);

        bool isSecondBestActive=myStr.optBaseIndex!=nodeNotActive;
//        if(isSecondBestActive){
//            std::cout<<"second best base "<<baseIDs.at(myStr.optBaseIndex)<<", lifted"<<std::endl;
//            for(auto it=secondBest.begin();it!=secondBest.end();it++){
//                std::cout<<*it<<",";
//            }
//            std::cout<<std::endl;
//        }
//        else{
//            std::cout<<"second best not active "<<std::endl;
//        }

        auto sbIt=secondBest.begin();
        listIt=isNotZeroInOpt.erase(listIt);

        //Comparing the list of optimal vertices with the second best solution and changing
        //isNotZeroInOpt and isOneInOpt accordingly
        while(listIt!=isNotZeroInOpt.end()&&sbIt!=secondBest.end()){

            if(*sbIt==*listIt){

                isOneInOpt.insert(*sbIt);
                sbIt++;
                listIt++;
            }
            else if(reachable(*sbIt,*listIt)){

                isOneInOpt.insert(*sbIt);
                //liftedMessages[*sbIt]=0;
                sbIt++;
            }
            else if(reachable(*listIt,*sbIt)){

                listIt=isNotZeroInOpt.erase(listIt);
            }
            else{

                listIt=isNotZeroInOpt.erase(listIt);
                isOneInOpt.insert(*sbIt);
                //liftedMessages[*sbIt]=0;
                sbIt++;
            }


        }

        isNotZeroInOpt.erase(listIt,isNotZeroInOpt.end());
        while(sbIt!=secondBest.end()){
            isOneInOpt.insert(*sbIt);
            //liftedMessages[*sbIt]=0;
            sbIt++;
        }


        listIt=isNotZeroInOpt.begin();

        double delta=currentOptValue-newOpt;

       // std::cout<<"delta "<<delta<<std::endl;

        size_t orderToClose=liftedIDToOrder.at(vertexToClose);
        localLiftedCosts[orderToClose]-=delta;
        liftedMessages[vertexToClose]=delta;
        minMarginalsImproving+=delta;
        currentOptValue=newOpt;


    }


    //Compute topDown structure and optimal solution after the change of lifted costs
//    myStr.topDownValuesStructure.fillWith(0);  //TODO implement this in reset in StrForUpdate
//    myStr.topDownVertexIDStructure.fillWith(getVertexToReach());
    fillWithValue<double>(ldpInstance.sncTDStructure,minVertex,maxVertex+1,0);
    fillWithValue<size_t>(ldpInstance.sncNeighborStructure,minVertex,maxVertex+1,getVertexToReach());
    topDownUpdate(myStr);
    ShiftedVector<char> verticesInScope(minVertex,maxVertex,ldpInstance.sncClosedVertices);


    assert(currentOptValue+minMarginalsImproving-origOptValue>-eps); //lower bound not decreasing
    assert(std::abs(currentOptValue-myStr.optValue)<eps);  //found optimal value as expected after changes of lifted costs

    //Structures for bottomUpUpdate
    //ShiftedVector<double> bottomUpValuesStructure(nodeID,mostDistantNeighborID,std::numeric_limits<double>::max());
    //ShiftedVector<char> closedVertices(nodeID,mostDistantNeighborID,false);

    fillWithValue<double>(ldpInstance.sncBUStructure,minVertex,maxVertex+1,std::numeric_limits<double>::max());
    fillWithValue<char>(ldpInstance.sncClosedVertices,minVertex,maxVertex+1,0);

    //The bottom up value for optimal vertices is known. It is obtained by subtracting the top down value of their descendants from the currentOptValue
    //Note that vertices closed in this for cycle will not have valid bottomUpVertexIDStructure entries
    for(size_t optVertex:isOneInOpt){
        assert(optVertex<ldpInstance.getNumberOfVertices());
        //size_t bestDesc=myStr.topDownVertexIDStructure[optVertex];
        size_t bestDesc=ldpInstance.sncNeighborStructure[optVertex];
        double toSubtract=0;
        //if(bestDesc!=getVertexToReach()) toSubtract=myStr.topDownValuesStructure[bestDesc];
        if(bestDesc!=getVertexToReach()) toSubtract=ldpInstance.sncTDStructure[bestDesc];
        ldpInstance.sncBUStructure[optVertex]=currentOptValue-toSubtract;
        ldpInstance.sncClosedVertices[optVertex]=1;
	}


    //Precompute some values for endpoints of base edges
   // ShiftedVector<size_t> bottomUpVertexIDStructure(nodeID,mostDistantNeighborID,getVertexToReach());
     //ShiftedVector<size_t> bottomUpVertexIDStructure(nodeID,mostDistantNeighborID,getVertexToReach());
     fillWithValue(ldpInstance.sncBUNeighborStructure,minVertex,maxVertex+1,getVertexToReach());
    for (int i = 0; i < baseIDs.size(); ++i) {
        if(baseIDs.at(i)==getVertexToReach()) continue;
        if(ldpInstance.sncClosedVertices[baseIDs.at(i)]) continue;
        ldpInstance.sncBUNeighborStructure[baseIDs.at(i)]=nodeID;
        ldpInstance.sncBUStructure[baseIDs.at(i)]=localBaseCosts.at(i);
    }

    //Compute min marginals for non-optimal nodes by finding best paths from the central node to these nodes
    //Traversing in bottom up order ensures that already computed values of bottom up structures remain valid
	for(size_t vertexID:liftedIDs){
        if(!ldpInstance.sncClosedVertices[vertexID]){
            std::unordered_map<size_t,double> newMessages=bottomUpUpdate(myStr,vertexID,verticesInScope,false);
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


    //Storing values from messages to output vector
	std::vector<double> messagesToOutput=std::vector<double>(liftedCosts.size());
	for (int i = 0; i < messagesToOutput.size(); ++i) {
		messagesToOutput[i]=liftedMessages[liftedIDs[i]];
	}

    if(debug()){
        StrForTopDownUpdate myStr2(localBaseCosts,localLiftedCosts,ldpInstance,nodeID,mostDistantNeighborID,getVertexToReach());
        topDownUpdate(myStr2);
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
       // if(debug()) std::cout<<"repam left "<<l.nodeID<<": "<<l.getLiftedID(right_node)<<":"<<msg<<std::endl;
		assert(msg_dim == 0);
        l.updateEdgeCost(msg,right_node,true);
	}

	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
	{
       // if(debug()) std::cout<<"repam right "<<r.nodeID<<": "<<r.getLiftedID(left_node)<<":"<<msg<<std::endl;
		assert(msg_dim == 0);
        r.updateEdgeCost(msg,left_node,true);
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
	{
        if(debug()) std::cout<<"message to left "<<r.nodeID<<": "<<r.getLiftedID(left_node)<<std::endl;
        const double delta = r.getOneLiftedMinMarginal(left_node);
		msg[0] -= omega * delta;
    //    if(debug()) std::cout<<"sent "<<r.nodeID<<": "<<r.getLiftedID(left_node)<<std::endl;
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_right(const SINGLE_NODE_CUT_FACTOR& l, MSG& msg, const double omega)
	{
        if(debug()) std::cout<<"message to right "<<l.nodeID<<": "<<l.getLiftedID(right_node)<<std::endl;
        const double delta = l.getOneLiftedMinMarginal(right_node);
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

//        std::vector<double> testMinMarginals=r.getAllLiftedMinMarginals();
//        for (int i = 0; i < testMinMarginals.size(); ++i) {
//            if(std::abs(testMinMarginals[i]-(1-omega)*msg_vec.at(i))>=eps){
//                std::cout<<"test min marginal "<<i<<": "<<testMinMarginals[i]<<", expected "<<(1-omega)*msg_vec.at(i)<<", omega "<<omega<<", orig min marginal "<<msg_vec.at(i)<<std::endl;
//                for (int j = 0; j < msg_vec.size(); ++j) {
//                    std::cout<<r.getLiftedID(j)<<": "<<testMinMarginals.at(j)<<", "<<msg_vec.at(j)<<std::endl;
//                }

//            }
//            assert(std::abs(testMinMarginals[i]-(1-omega)*msg_vec.at(i))<eps);
//        }

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

//        std::vector<double> testMinMarginals=l.getAllLiftedMinMarginals();
//        for (int i = 0; i < testMinMarginals.size(); ++i) {

//            if(std::abs(testMinMarginals[i]-(1-omega)*msg_vec.at(i))>=eps){
//                std::cout<<"test min marginal "<<i<<": "<<testMinMarginals[i]<<", expected "<<(1-omega)*msg_vec.at(i)<<", omega "<<omega<<", orig min marginal "<<msg_vec.at(i)<<std::endl;
//                for (int j = 0; j < msg_vec.size(); ++j) {
//                    std::cout<<l.getLiftedID(j)<<": "<<testMinMarginals.at(j)<<", "<<msg_vec.at(j)<<std::endl;
//                }

//            }
//            assert(std::abs(testMinMarginals[i]-(1-omega)*msg_vec.at(i))<eps);
//        }

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
