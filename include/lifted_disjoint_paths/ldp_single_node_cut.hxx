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



struct StrForTopDownUpdate{

    //Structure storing results of topDownUpdate
    StrForTopDownUpdate(const std::vector<double>& bCosts,const std::vector<double>& lCosts):
	baseCosts(bCosts),
    liftedCosts(lCosts),
    solutionCosts(bCosts.size()+1),
	optValue(0),
    //nodeID(centralNodeID),
    optBaseIndex(bCosts.size())
    {}

   // const size_t nodeID;
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

    void initNodeCost(double fractionNode);

	const std::vector<double>& getLiftedCosts() const {
		return liftedCosts;
    }

	const std::vector<double>& getBaseCosts() const {
		return baseCosts;
	}

    const double& getNodeCost()const{
        return nodeCost;
    }

    //Methods for getting and setting primal solution
    double EvaluatePrimal() const;

	void init_primal(){
		primalBase_=nodeNotActive;
	}

	void setBaseEdgeActive(size_t index);
	void setNoBaseEdgeActive();

    void setPrimalLifted(std::vector<size_t> &verticesOfActiveEdges);

	bool isNodeActive() const
	{
		return primalBase_!=nodeNotActive;
	}

    bool isNodeOutFlow() const{
        return isOutFlow;
    }

	size_t getPrimalBaseIndex() const {
		return primalBase_;
	}


	size_t getPrimalBaseVertexID() const {
        assert(primalBase_<baseIDs.size());
        return baseIDs.at(primalBase_);
	}

    const std::vector<size_t>& getPrimalLiftedIndices() const {
		return primalLifted_;
	}

//	const bool isActiveInPrimalLifted(size_t vertex) const {
//		return primalLifted_.count(vertex)>0;
//		//return primalLifted_.at(vertex);
//	}


    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(baseCosts, liftedCosts); }

    auto export_variables() { return std::tie(baseCosts, liftedCosts); }

    //Updating costs
    void updateEdgeCost(const double value,const size_t vertexIndex,bool isLifted);
	void updateNodeCost(const double value);


    //Computing lower bound and min marginals
    double LowerBound() const;
    double LowerBound(bool canChange) const;

    double getOneBaseEdgeMinMarginal(size_t vertex, const std::vector<double>*pBaseCosts, const std::vector<double>*pLiftedCosts) const;
    double getOneLiftedMinMarginal(size_t indexOfLiftedEdge, const std::vector<double>*pBaseCosts, const std::vector<double>*pLiftedCosts)const;
    double getNodeMinMarginal()const;

    std::vector<double> getAllBaseMinMarginals(const std::vector<double> *pLocalBaseCosts, const std::vector<double> *pLocalLiftedCosts)const;
    std::vector<double> getAllBaseMinMarginals()const;
    std::vector<double> getAllBaseMinMarginalsForMCF() const;
    std::vector<double> getAllLiftedMinMarginals(const std::vector<double> *pLocalBaseCosts=nullptr,const std::vector<double> *pLocalLiftedCosts=nullptr) const;


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

    const double& getPrimalLiftedCost()const{
        return primalLiftedCost;
    }
    const double& getPrimalBaseCost()const{
        return primalBaseCost;
    }


         const std::size_t nodeID;


private:
    //Methods used in computing lower bound or min marginals
    void topDownUpdate(StrForTopDownUpdate& myStr, const size_t vertexIDToIgnore) const;
    void topDownUpdate(StrForTopDownUpdate& myStr) const;
    void bottomUpUpdate(const StrForTopDownUpdate& myStr, const size_t vertex)const;
    void updateOptimal() const;
    void initTraverseOrder();

    //Obtain IDs of vertices of lifted edges that are part of a found optimal solution
    std::list<size_t> getOptLiftedFromIndexStr(const StrForTopDownUpdate& myStr)const;


    //Methods for exploring graph structures
    bool isInGivenInterval(const size_t nodeIndex,const size_t boundaryIndex) const {
        assert(nodeIndex < ldpInstance.getNumberOfVertices());
        if(isOutFlow){
            return nodeIndex<=boundaryIndex;
        }
        else{
            return nodeIndex>=boundaryIndex;
        }
    }

    const std::pair<size_t,double>* neighborsBegin(const size_t& nodeIndex)const{
        if(isOutFlow){
            return ldpBaseGraph.forwardNeighborsBegin(nodeIndex);
        }
        else {
            return ldpBaseGraph.backwardNeighborsBegin(nodeIndex);
        }
    }

    const std::pair<size_t,double>* neighborsEnd(const size_t& nodeIndex)const{
        if(isOutFlow){
            return ldpBaseGraph.forwardNeighborsEnd(nodeIndex);
        }
        else {
            return ldpBaseGraph.backwardNeighborsEnd(nodeIndex);
        }
    }

    const std::pair<size_t,double>* neighborsRevBegin(const size_t& nodeIndex)const{
        if(!isOutFlow){
            return ldpBaseGraph.forwardNeighborsBegin(nodeIndex);
        }
        else {
            return ldpBaseGraph.backwardNeighborsBegin(nodeIndex);
        }
    }

    const std::pair<size_t,double>* neighborsRevEnd(const size_t& nodeIndex)const{
        if(!isOutFlow){
            return ldpBaseGraph.forwardNeighborsEnd(nodeIndex);
        }
        else {
            return ldpBaseGraph.backwardNeighborsEnd(nodeIndex);
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
     const LDP_INSTANCE& ldpInstance;
     const LdpDirectedGraph& ldpBaseGraph;
     const LdpDirectedGraph& ldpLiftedGraph;

     //costs
     std::vector<double> baseCosts;
     std::vector<double> liftedCosts;
     double nodeCost;

     //Variables storing primal solutions
     size_t primalBase_;
     std::vector<size_t> primalLifted_;

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

     std::vector<size_t> traverseOrder;

     mutable double primalBaseCost;
     mutable double primalLiftedCost;



};


template<class LDP_INSTANCE>
inline  ldp_single_node_cut_factor<LDP_INSTANCE>::ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut):
nodeID(nID),
ldpInstance(ldpInst),
isOutFlow(isOut),
ldpBaseGraph(ldpInst.getMyGraph()),
ldpLiftedGraph(ldpInst.getMyGraphLifted())

{


    mostDistantNeighborID=nodeID;
    initBaseCosts(0);
    initLiftedCosts(0);
    initNodeCost(0);

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

	optValue=0;

    nodeCost=0;
    optValueUpToDate=false;
    solutionCostsUpToDate=false;


    initTraverseOrder();

}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::initTraverseOrder() {


        fillWithValue<char>(ldpInstance.sncClosedVertices,minVertex,maxVertex+1,0);



    std::stack<size_t> nodeStack;
    nodeStack.push(nodeID);


    traverseOrder.clear();
    while(!nodeStack.empty()){
        size_t currentNode=nodeStack.top();

        assert(currentNode<ldpInstance.getNumberOfVertices());
        if(ldpInstance.sncClosedVertices[currentNode]){
            nodeStack.pop();
        }
        else{

            bool descClosed=true;
             //Traverse all base neighbors in direction out from central node
            //Not closed vertices store to the stack.
            const std::pair<size_t,double>* vertexIt=neighborsBegin(currentNode);
            const std::pair<size_t,double>* end=neighborsEnd(currentNode);


            for (;vertexIt!=end;vertexIt++) {


               size_t desc=vertexIt->first;

              if(desc==getVertexToReach()) continue;

                if(isInGivenInterval(desc,mostDistantNeighborID)){

                    if(!ldpInstance.sncClosedVertices[desc]){  //descendant closed
                        nodeStack.push(desc);
                        descClosed=false;
                    }
                }


            }
            if(descClosed){
                traverseOrder.push_back(currentNode);
                assert(currentNode<ldpInstance.getNumberOfVertices());
                ldpInstance.sncClosedVertices[currentNode]=1;
                nodeStack.pop();

            }
        }

    }

}

template<class LDP_INSTANCE>
inline std::list<size_t> ldp_single_node_cut_factor<LDP_INSTANCE>::getOptLiftedFromIndexStr(const StrForTopDownUpdate& myStr) const{

	std::list<size_t> optLifted;
    double optValueComputed=0;
    if(myStr.optBaseIndex!=nodeNotActive){
        optValueComputed=myStr.baseCosts.at(myStr.optBaseIndex);
        optValueComputed+=nodeCost;

        size_t vertexInOptimalPath=baseIDs.at(myStr.optBaseIndex);

        while(vertexInOptimalPath!=getVertexToReach()){

            if(isLiftedVertex(vertexInOptimalPath)){

                optLifted.push_back(vertexInOptimalPath);
                if(debug()){
                    double toAdd=myStr.liftedCosts.at(liftedIDToOrder.at(vertexInOptimalPath));
                    optValueComputed+=toAdd;
                }

            }

            assert(vertexInOptimalPath<ldpInstance.getNumberOfVertices());
            vertexInOptimalPath=ldpInstance.sncNeighborStructure[vertexInOptimalPath];
		}
        if(debug()) assert(std::abs(optValueComputed-myStr.optValue)<eps);
    }

    return optLifted;

}





template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::setPrimalLifted(std::vector<size_t> &verticesOfActiveEdges) {
	primalLifted_=verticesOfActiveEdges;
}





template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::EvaluatePrimal() const{
    primalBaseCost=0;
    primalLiftedCost=0;
    if(primalBase_==nodeNotActive){

        return 0;
    }
    else{
       // std::cout<<"node cost "<<nodeID<<": "<<nodeCost<<std::endl;
        double value=nodeCost;
       // std::cout<<"node cost "<<nodeID<<": "<<nodeCost<<std::endl;
        value+=baseCosts.at(primalBase_);
        primalBaseCost=value;

        for(size_t node:primalLifted_){
            value+=liftedCosts.at(node);
            primalLiftedCost+=liftedCosts.at(node);
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
            assert(solutionCosts.size()>0);
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
            StrForTopDownUpdate myStr(baseCosts,liftedCosts);
            topDownUpdate(myStr);
            optValue=myStr.optValue;
            solutionCosts=myStr.solutionCosts;
            optBaseIndex=myStr.optBaseIndex;
            solutionCostsUpToDate=true;


        }
        optValueUpToDate=true;
    }



}

template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getNodeMinMarginal()const{

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
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getOneBaseEdgeMinMarginal(const size_t index, const std::vector<double>*pBaseCosts, const std::vector<double>*pLiftedCosts)const{
	assert(index<baseCosts.size());
    assert(optBaseIndex<solutionCosts.size());

    StrForTopDownUpdate strForUpdateValues(*pBaseCosts,*pLiftedCosts);
    topDownUpdate(strForUpdateValues);


    if(strForUpdateValues.optBaseIndex!=index){
        //return strForUpdateValues.solutionCosts[index]-strForUpdateValues.solutionCosts[strForUpdateValues.optBaseIndex];
        return strForUpdateValues.solutionCosts[index]-strForUpdateValues.optValue;
    }
    else{
        double secondBest=std::numeric_limits<double>::max();
        for (int i = 0; i < strForUpdateValues.solutionCosts.size(); ++i) {
            if(i==strForUpdateValues.optBaseIndex) continue;
            if(strForUpdateValues.solutionCosts[i]<secondBest){
                secondBest=strForUpdateValues.solutionCosts[i];
            }
        }
        return strForUpdateValues.solutionCosts[index]-secondBest;
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
inline std::vector<double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllBaseMinMarginals(const std::vector<double>* pLocalBaseCosts,const std::vector<double>* pLocalLiftedCosts) const{


    StrForTopDownUpdate str(*pLocalBaseCosts,*pLocalLiftedCosts);
    topDownUpdate(str);

    std::vector<double> minMarginals(pLocalBaseCosts->size());
    if(str.optBaseIndex==nodeNotActive){
        for (int i = 0; i < str.solutionCosts.size()-1; ++i) {
            minMarginals[i]=str.solutionCosts[i];

        }
    }
    else{
        double secondBest=std::numeric_limits<double>::infinity();
        double optValue=str.solutionCosts.at(str.optBaseIndex);
        for (int i = 0; i < str.solutionCosts.size(); ++i) {
            if(i==str.optBaseIndex) continue;
            if(str.solutionCosts[i]<secondBest){
                secondBest=str.solutionCosts[i];
            }

        }
        for (int i = 0; i < str.solutionCosts.size()-1; ++i) {
            minMarginals[i]=str.solutionCosts[i]-secondBest;

        }
    }
    return minMarginals;
}

template<class LDP_INSTANCE>
inline std::vector<double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllBaseMinMarginals() const{
    return getAllBaseMinMarginals(&baseCosts,&liftedCosts);
//	updateOptimal();

//    std::vector<double> minMarginals(baseCosts.size());
//    if(optBaseIndex==nodeNotActive){
//		for (int i = 0; i < solutionCosts.size()-1; ++i) {
//			minMarginals[i]=solutionCosts[i];

//		}
//	}
//	else{
//		double secondBest=std::numeric_limits<double>::infinity();
//        double optValue=solutionCosts.at(optBaseIndex);
//		for (int i = 0; i < solutionCosts.size(); ++i) {
//            if(i==optBaseIndex) continue;
//			if(solutionCosts[i]<secondBest){
//				secondBest=solutionCosts[i];
//			}

//		}
//		for (int i = 0; i < solutionCosts.size()-1; ++i) {
//			minMarginals[i]=solutionCosts[i]-secondBest;

//		}
//    }
//	return minMarginals;
}

template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateNodeCost(const double value){
     // std::cout<<"update node cost "<<nodeID<<", out "<<isOutFlow<<": "<<value<<std::endl;
	nodeCost+=value;
	for (int i = 0; i < solutionCosts.size()-1; ++i) {
		(solutionCosts[i]+=value);
	}
    optValueUpToDate=false;

}

template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateEdgeCost(const double value,const size_t vertexIndex,bool isLifted){//Only cost change
  //  if(debug()) std::cout<<"updating cost "<<value<<std::endl;
	if(!isLifted){ //update in base edge
		assert(vertexIndex<baseCosts.size());
       // std::cout<<"update edge, node id "<<nodeID<<", "<<baseIDs[vertexIndex]<<": "<<value<<std::endl;
		baseCosts[vertexIndex]+=value;
		solutionCosts[vertexIndex]+=value;
        optValueUpToDate=false;
      //  std::cout<<"update edge cost "<<nodeID<<" "<<baseIDs[vertexIndex]<<": "<<value<<std::endl;
	}
	else{ //update in lifted edge

		assert(vertexIndex<liftedCosts.size());
		liftedCosts[vertexIndex]+=value;

        optValueUpToDate=false;
        solutionCostsUpToDate=false;


	}
    //LowerBound();
}



template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::initNodeCost(double fractionNode){
    nodeCost=fractionNode*ldpInstance.getVertexScore(nodeID);

}

template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::initLiftedCosts(double fractionLifted){
	liftedCosts=std::vector<double>();
	liftedIDs=std::vector<size_t>();
	liftedIDToOrder.clear();
    const LdpDirectedGraph& myLiftedGraph=ldpInstance.getMyGraphLifted();

    if(isOutFlow){
        const std::pair<size_t,double>* edgeIt=myLiftedGraph.forwardNeighborsBegin(nodeID);

        size_t counter=0;
        for (;edgeIt!=myLiftedGraph.forwardNeighborsEnd(nodeID);edgeIt++) {
            const size_t& node=edgeIt->first;
            const double& cost=edgeIt->second;

            liftedCosts.push_back((cost)*fractionLifted);
            liftedIDs.push_back(node);
            liftedIDToOrder[node]=counter;
            if(node!=ldpInstance.getTerminalNode()){
                mostDistantNeighborID=std::max(mostDistantNeighborID,node);
            }
            counter++;
        }
    }
    else{
        const std::pair<size_t,double>* edgeIt=myLiftedGraph.backwardNeighborsBegin(nodeID);
        size_t counter=0;
        for (;edgeIt!=myLiftedGraph.backwardNeighborsEnd(nodeID);edgeIt++) {
            const size_t& node=edgeIt->first;
            const double& cost=edgeIt->second;

            liftedCosts.push_back((cost)*fractionLifted);
            liftedIDs.push_back(node);
            liftedIDToOrder[node]=counter;
            if(node!=ldpInstance.getSourceNode()){
                mostDistantNeighborID=std::min(mostDistantNeighborID,node);
            }
            counter++;
        }
    }


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

        const std::pair<size_t,double>* edgeIt=myBaseGraph.forwardNeighborsBegin(nodeID);
        size_t counter=0;
        for (;edgeIt!=myBaseGraph.forwardNeighborsEnd(nodeID);edgeIt++) {
            const size_t& node=edgeIt->first;
            const double& cost=edgeIt->second;

            baseCosts.push_back((cost)*fractionBase);

            baseIDs.push_back(node);
            baseIDToIndex[node]=counter;
            if(node!=ldpInstance.getTerminalNode()){
                mostDistantNeighborID=std::max(mostDistantNeighborID,node);
            }
            counter++;
        }
    }
    else{
        const std::pair<size_t,double>* edgeIt=myBaseGraph.backwardNeighborsBegin(nodeID);
        size_t counter=0;
        for (;edgeIt!=myBaseGraph.backwardNeighborsEnd(nodeID); edgeIt++) {
            const size_t& node=edgeIt->first;
            const double& cost=edgeIt->second;

            baseCosts.push_back((cost)*fractionBase);
            //std::cout<<"snc base in "<<nodeID<<" "<<node<<": "<<cost<<std::endl;
            baseIDs.push_back(node);
            baseIDToIndex[node]=counter;
            if(node!=ldpInstance.getSourceNode()){
                mostDistantNeighborID=std::min(mostDistantNeighborID,node);
            }
            counter++;
        }
    }




    solutionCosts=std::vector<double>(baseCosts.size()+1);

    optValueUpToDate=false;
    solutionCostsUpToDate=false;

}



template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::LowerBound(bool canChange) const{//TODO store info about how valuesStructures changed. At least max time layer of changed lifted edge
    //if(canChange||optValueUpToDate){
    if(canChange){
        updateOptimal();
        if(debug()){
            double controlValue=0;
            if(optBaseIndex!=nodeNotActive){
                controlValue+=nodeCost;
                controlValue+=baseCosts.at(optBaseIndex);

            }
        }

        //std::cout<<"snc lower bound "<<optValue<<std::endl;
        return optValue;
    }
    else{
        StrForTopDownUpdate myStr(baseCosts,liftedCosts);
        topDownUpdate(myStr);
        return myStr.optValue;
    }

}


template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::LowerBound() const{//TODO store info about how valuesStructures changed. At least max time layer of changed lifted edge
    return LowerBound(false);

}





//TODO make this separate from standard update
template<class LDP_INSTANCE>
void ldp_single_node_cut_factor<LDP_INSTANCE>::topDownUpdate(StrForTopDownUpdate& myStr) const{
     topDownUpdate(myStr,getVertexToReach());
}


template<class LDP_INSTANCE>
void ldp_single_node_cut_factor<LDP_INSTANCE>::topDownUpdate(StrForTopDownUpdate& myStr,const size_t vertexIDToIgnore) const{


    bool vertexToIgnoreSet=false;
    size_t lastVertex=mostDistantNeighborID;
    if(vertexIDToIgnore!=getVertexToReach()){
        vertexToIgnoreSet=true;
        lastVertex=vertexIDToIgnore;
        fillWithValue<double>(ldpInstance.sncTDStructure,std::min(nodeID,vertexIDToIgnore),std::max(nodeID,vertexIDToIgnore)+1,0);

    }
    else{
        for(size_t v:traverseOrder){
            ldpInstance.sncTDStructure[v]=0;
            ldpInstance.sncNeighborStructure[v]=getVertexToReach();

        }
    }



    //Store all lifted costs to top down values structure
    for (int i = 0; i < liftedIDs.size(); ++i) {
        if(!vertexToIgnoreSet||isInGivenInterval(liftedIDs.at(i),lastVertex)){
            ldpInstance.sncTDStructure[liftedIDs.at(i)]=myStr.liftedCosts.at(i);
        }
    }




    for (size_t i=0; i < traverseOrder.size()-1; ++i) {

        size_t currentNode=traverseOrder[i];

        if(currentNode==vertexIDToIgnore) continue;
        if(!isInGivenInterval(currentNode,lastVertex)) continue;
        double bestDescValue=0;
        size_t bestDescVertexID=getVertexToReach();

         //Search for best descendant
        const std::pair<size_t,double>* vertexIt=neighborsBegin(currentNode);
        const std::pair<size_t,double>* end=neighborsEnd(currentNode);


        for (;vertexIt!=end;vertexIt++) {

            size_t desc=vertexIt->first;

            if(desc==vertexIDToIgnore||desc==getVertexToReach()) continue;

            if(isInGivenInterval(desc,mostDistantNeighborID)){

                double value=ldpInstance.sncTDStructure[desc];
                if(bestDescValue>value){
                    bestDescValue=value;
                    bestDescVertexID=desc;
                }
            }
        }



        ldpInstance.sncTDStructure[currentNode]+=bestDescValue;
        ldpInstance.sncNeighborStructure[currentNode]=bestDescVertexID;

    }


    //all nodes closed, compute solution values
    double bestSolutionValue=0;
    size_t bestSolutionIndex=nodeNotActive;

    myStr.solutionCosts[nodeNotActive]=0;
    for (size_t i = 0; i < myStr.baseCosts.size(); ++i) {
        double baseCost=myStr.baseCosts[i];
        if(!vertexToIgnoreSet||baseIDs[i]!=vertexIDToIgnore){

            double valueToAdd=0;
            if(baseIDs[i]!=getVertexToReach()){
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







template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getOneLiftedMinMarginal(size_t indexOfLiftedEdge, const std::vector<double> *pBaseCosts, const std::vector<double> *pLiftedCosts)const{
    assert(indexOfLiftedEdge<liftedCosts.size());

    // std::cout<<"one lifted min marginal in snc"<<std::endl;

    const std::vector<double>&localBaseCosts=*pBaseCosts;
    const std::vector<double>&localLiftedCosts=*pLiftedCosts;

    StrForTopDownUpdate strForUpdateValues(localBaseCosts,localLiftedCosts);
    topDownUpdate(strForUpdateValues);
    double origOptValue=strForUpdateValues.optValue;


    std::list<size_t> optimalSolutionLifted=getOptLiftedFromIndexStr(strForUpdateValues);

	bool isOptimal=false;
	for(size_t optVertex:optimalSolutionLifted){
        if(optVertex==liftedIDs.at(indexOfLiftedEdge)){
			isOptimal=true;

			break;
		}
	}


	if(isOptimal){

        //If is optimal, run topDownUpdate with ignoring this vertex
        topDownUpdate(strForUpdateValues,liftedIDs.at(indexOfLiftedEdge));
        double restrictedOptValue=strForUpdateValues.optValue;

        double valueToReturn=origOptValue-restrictedOptValue;
        assert(valueToReturn<eps);
        if(valueToReturn>=eps){
            throw std::runtime_error("wrong lifted message opt");
        }
		return valueToReturn;



	}
	else{

        //If it is not optimal, find the best possible solution containing this vertex active via bottomUpUpdate
       for(size_t v:traverseOrder){
            ldpInstance.sncBUStructure[v]=std::numeric_limits<double>::max();
            ldpInstance.sncClosedVertices[v]=0;
            ldpInstance.sncBUNeighborStructure[v]=getVertexToReach();
        }


        for (int i = 0; i < baseIDs.size(); ++i) {
            if(baseIDs.at(i)==getVertexToReach()) continue;
            ldpInstance.sncBUNeighborStructure[baseIDs.at(i)]=nodeID;
            ldpInstance.sncBUStructure[baseIDs.at(i)]=localBaseCosts.at(i);
        }

//        ShiftedVector<char> verticesInScope(minVertex,maxVertex);
//        for(size_t i=0;i<traverseOrder.size();i++){
//            verticesInScope[traverseOrder[i]]=1;
//        }
        bottomUpUpdate(strForUpdateValues,liftedIDs[indexOfLiftedEdge]);

        //auto it =message.begin();
        double messValue=ldpInstance.sncLiftedMessages[liftedIDs[indexOfLiftedEdge]];

        assert(messValue>-eps);
        if(messValue<=-eps){
            throw std::runtime_error("wrong lifted message");
        }

        return messValue;
	}

}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::bottomUpUpdate(const StrForTopDownUpdate& myStr,const size_t vertex)const{


    bool onlyOne=vertex!=nodeID;
    for(size_t i=0;i<traverseOrder.size();i++){
        ldpInstance.sncVerticesInScope[traverseOrder[i]]=1;
    }


    for(size_t i=traverseOrder.size();i>=1;i--){
        size_t currentVertex=traverseOrder[i-1];
        if(ldpInstance.sncClosedVertices[currentVertex]) continue;
        size_t bestIndex=getVertexToReach();
        double bestValue=std::numeric_limits<double>::max();

        if(ldpInstance.sncBUNeighborStructure[currentVertex]==nodeID){  //This holds for endpoints of base edges
            bestValue=ldpInstance.sncBUStructure[currentVertex]+nodeCost;
            bestIndex=nodeID;
        }
        const std::pair<size_t,double>* vertexIt=neighborsRevBegin(currentVertex);
        const std::pair<size_t,double>* end=neighborsRevEnd(currentVertex);
        for (;vertexIt!=end;vertexIt++) {

            size_t pred=vertexIt->first;
            assert(pred<ldpInstance.getNumberOfVertices());

            bool newConstraint=(pred==nodeID||!ldpInstance.sncVerticesInScope[pred]);
            if(newConstraint) continue;
            assert(ldpInstance.sncClosedVertices[pred]>0);
            double value=ldpInstance.sncBUStructure[pred];
            if(value<bestValue){
                bestValue=value;
                bestIndex=pred;  //TODO check
            }

        }

        ldpInstance.sncBUNeighborStructure[currentVertex]=bestIndex;
        if(isLiftedVertex(currentVertex)){
            if(onlyOne&&currentVertex!=vertex){  //For the case of getting one min marginal
                assert(bestValue!=std::numeric_limits<double>::max());
                bestValue+=myStr.liftedCosts.at(liftedIDToOrder.at(currentVertex));
            }
            else{
                //
                assert(bestValue!=std::numeric_limits<double>::max());
                double topDownValueOfDesc=0;

                size_t bestDesc=ldpInstance.sncNeighborStructure[currentVertex]; //Best neighbor in the direction from the central node
                if(bestDesc!=getVertexToReach()){
                    topDownValueOfDesc=ldpInstance.sncTDStructure[bestDesc];
                }

                double topDownCurrent=ldpInstance.sncTDStructure[currentVertex];  //Best top down value

                double restrictedOpt=topDownCurrent+bestValue;  //Best top down value of currentVertex plus bottom up value of best predecessor
                double delta=restrictedOpt-myStr.optValue;
                bestValue=myStr.optValue-topDownValueOfDesc;  //Compute node's bottom up value after changing its lifted cost by delta

                ldpInstance.sncLiftedMessages[currentVertex]=delta;

                if(onlyOne) break;

            }
        }

        ldpInstance.sncClosedVertices[currentVertex]=1;
        ldpInstance.sncBUStructure[currentVertex]=bestValue;

    }
    for(size_t i=0;i<traverseOrder.size();i++){
        ldpInstance.sncVerticesInScope[traverseOrder[i]]=0;
    }


}


template<class LDP_INSTANCE>
inline std::vector<double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllLiftedMinMarginals(const std::vector<double>* pLocalBaseCosts, const std::vector<double> *pLocalLiftedCosts) const{



    std::vector<double> localLiftedCosts;
    std::vector<double> localBaseCosts;
    if(pLocalBaseCosts==nullptr){
        localBaseCosts=baseCosts;
    }
    else{//If other value of base costs needs to be used. E.g. for getting min marginals of both base and lifted edges
        localBaseCosts=*pLocalBaseCosts;
    }


    if(pLocalLiftedCosts==nullptr){
        localLiftedCosts=liftedCosts;
    }
    else{
        localLiftedCosts=*pLocalLiftedCosts;
    }

    //First, compute optimal value
    StrForTopDownUpdate myStr(localBaseCosts,localLiftedCosts);
    topDownUpdate(myStr);
    double origOptValue=myStr.optValue;


    //All vertices that are not zero in any optimal solution
    std::list<size_t> isNotZeroInOpt=getOptLiftedFromIndexStr(myStr);
    //All vertices that are one in at least one of the optimal solutions
    std::unordered_set<size_t> isOneInOpt(isNotZeroInOpt.begin(),isNotZeroInOpt.end());



    double minMarginalsImproving=0;

    double currentOptValue=myStr.optValue;

    auto listIt=isNotZeroInOpt.begin();


    for(size_t v: liftedIDs){
        ldpInstance.sncLiftedMessages[v]=0;
    }

    //Obtaining min marginals for nodes that are active in all optimal solutions
    while(!isNotZeroInOpt.empty()){
        size_t vertexToClose=*listIt;


        //Obtaining best solution while ignoring vertexToClose
        topDownUpdate(myStr,vertexToClose);
        double newOpt=myStr.optValue;
        std::list<size_t> secondBest=getOptLiftedFromIndexStr(myStr);

        bool isSecondBestActive=myStr.optBaseIndex!=nodeNotActive;

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

        size_t orderToClose=liftedIDToOrder.at(vertexToClose);
        localLiftedCosts[orderToClose]-=delta;
        ldpInstance.sncLiftedMessages[vertexToClose]=delta;
        minMarginalsImproving+=delta;
        currentOptValue=newOpt;


    }


    //Compute topDown structure and optimal solution after the change of lifted costs
    topDownUpdate(myStr);
//    ShiftedVector<char> verticesInScope(minVertex,maxVertex);
//    for(size_t i=0;i<traverseOrder.size();i++){
//        verticesInScope[traverseOrder[i]]=1;
//    }



    assert(currentOptValue+minMarginalsImproving-origOptValue>-eps); //lower bound not decreasing
    assert(std::abs(currentOptValue-myStr.optValue)<eps);  //found optimal value as expected after changes of lifted costs

    //Structures for bottomUpUpdate
    for(size_t v:traverseOrder){
        ldpInstance.sncBUStructure[v]=std::numeric_limits<double>::max();
        ldpInstance.sncClosedVertices[v]=0;
        ldpInstance.sncBUNeighborStructure[v]=getVertexToReach();
    }

    //The bottom up value for optimal vertices is known. It is obtained by subtracting the top down value of their descendants from the currentOptValue
    //Note that vertices closed in this for cycle will not have valid bottomUpVertexIDStructure entries
    for(size_t optVertex:isOneInOpt){
        assert(optVertex<ldpInstance.getNumberOfVertices());
        size_t bestDesc=ldpInstance.sncNeighborStructure[optVertex];
        double toSubtract=0;
        if(bestDesc!=getVertexToReach()) toSubtract=ldpInstance.sncTDStructure[bestDesc];
        ldpInstance.sncBUStructure[optVertex]=currentOptValue-toSubtract;
        ldpInstance.sncClosedVertices[optVertex]=1;

	}


    //Precompute some values for endpoints of base edges
    for (int i = 0; i < baseIDs.size(); ++i) {
        if(baseIDs.at(i)==getVertexToReach()) continue;
        if(ldpInstance.sncClosedVertices[baseIDs.at(i)]) continue;
        ldpInstance.sncBUNeighborStructure[baseIDs.at(i)]=nodeID;
        ldpInstance.sncBUStructure[baseIDs.at(i)]=localBaseCosts.at(i);
    }

    //Compute min marginals for non-optimal nodes by finding best paths from the central node to these nodes
    //Traversing in bottom up order ensures that already computed values of bottom up structures remain valid
    bottomUpUpdate(myStr,nodeID);


    //Storing values from messages to output vector
    std::vector<double> messagesToOutput=std::vector<double>(liftedCosts.size());
	for (int i = 0; i < messagesToOutput.size(); ++i) {
        messagesToOutput[i]=ldpInstance.sncLiftedMessages[liftedIDs[i]];
	}

    if(debug()){
        StrForTopDownUpdate myStr2(localBaseCosts,localLiftedCosts);
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
       // if(debug()) l.LowerBound();
		assert(msg_dim == 0);
        l.updateEdgeCost(msg,right_node,true);
       // if(debug()) l.LowerBound();
	}

	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
	{
       // if(debug()) std::cout<<"repam right "<<r.nodeID<<": "<<r.getLiftedID(left_node)<<":"<<msg<<std::endl;
      // if(debug()) r.LowerBound();
        assert(msg_dim == 0);
        r.updateEdgeCost(msg,left_node,true);
      //  if(debug()) r.LowerBound();
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
	{
        const std::vector<double>& baseCosts=r.getBaseCosts();
        const std::vector<double>& liftedCosts=r.getLiftedCosts();
        //if(debug()) std::cout<<"message to left "<<r.nodeID<<": "<<r.getLiftedID(left_node)<<std::endl;
        const double delta = r.getOneLiftedMinMarginal(left_node,&baseCosts,&liftedCosts);
        if(debug()){
            std::vector<double> controlLiftedCosts=liftedCosts;
            controlLiftedCosts[left_node]-=delta;
            double controlDelta=r.getOneLiftedMinMarginal(left_node,&baseCosts,&controlLiftedCosts);
            if(abs(controlDelta)>eps){
                throw std::runtime_error("different min marginal after lower bound");
            }
        }
		msg[0] -= omega * delta;
    //    if(debug()) std::cout<<"sent "<<r.nodeID<<": "<<r.getLiftedID(left_node)<<std::endl;
	}

	template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
	void send_message_to_right(const SINGLE_NODE_CUT_FACTOR& l, MSG& msg, const double omega)
	{
        const std::vector<double>& baseCosts=l.getBaseCosts();
        const std::vector<double>& liftedCosts=l.getLiftedCosts();
       // if(debug()) std::cout<<"message to right "<<l.nodeID<<": "<<l.getLiftedID(right_node)<<std::endl;
        const double delta = l.getOneLiftedMinMarginal(right_node,&baseCosts,&liftedCosts);
        if(debug()){
            std::vector<double> controlLiftedCosts=liftedCosts;
            controlLiftedCosts[right_node]-=delta;
            double controlDelta=l.getOneLiftedMinMarginal(right_node,&baseCosts,&controlLiftedCosts);
            if(abs(controlDelta)>eps){
                throw std::runtime_error("different min marginal after lower bound");
            }
        }
		msg[0] -= omega * delta;
	}



    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG_ARRAY>
    static void SendMessagesToLeft(const SINGLE_NODE_CUT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
    {
        const std::vector<double> msg_vec = r.getAllLiftedMinMarginals();

        if(debug()){
            std::vector<double> liftedCosts=r.getLiftedCosts();
            std::vector<double> baseCosts=r.getBaseCosts();

            for (int i = 0; i < liftedCosts.size(); ++i) {
                 liftedCosts[i]-=msg_vec[i];
            }

            const std::vector<double> controlMessages=r.getAllLiftedMinMarginals(&baseCosts,&liftedCosts);

            for (int i = 0; i < controlMessages.size(); ++i) {


                if(abs(controlMessages.at(i))>eps){
                    std::cout<<"WRONG control message "<<controlMessages.at(i)<<", orig message "<<msg_vec.at(i)<<std::endl;
                    throw std::runtime_error("wrong lifted min marginal in snc lifted message");
                }
            }

        }

        for(auto it=msg_begin; it!=msg_end; ++it)
        {

            auto& msg = (*it).GetMessageOp();
            const size_t left_node = msg.left_node;
            const size_t right_node = msg.right_node;
            double delta=omega * msg_vec.at(left_node);

            (*it)[0] -= delta;
        }





    }

    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG_ARRAY>
    static void SendMessagesToRight(const SINGLE_NODE_CUT_FACTOR& l, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
    {
        //if(debug()) std::cout<<"running get all lifted marginals to right "<<l.nodeID<<std::endl;

        const std::vector<double> msg_vec = l.getAllLiftedMinMarginals();

        if(debug()){
            std::vector<double> liftedCosts=l.getLiftedCosts();
            std::vector<double> baseCosts=l.getBaseCosts();

            for (int i = 0; i < liftedCosts.size(); ++i) {

                 liftedCosts[i]-=msg_vec[i];

             }

            const std::vector<double> controlMessages=l.getAllLiftedMinMarginals(&baseCosts,&liftedCosts);

            for (int i = 0; i < controlMessages.size(); ++i) {

                if(abs(controlMessages.at(i))>eps){
                    std::cout<<"WRONG control message "<<controlMessages.at(i)<<", orig message "<<msg_vec.at(i)<<std::endl;
                    throw std::runtime_error("wrong lifted min marginal in snc lifted message");
                }
            }

        }


        for(auto it=msg_begin; it!=msg_end; ++it)
        {

            auto& msg = (*it).GetMessageOp();
            const size_t left_node = msg.left_node;
            const size_t right_node = msg.right_node;

            double delta= omega * msg_vec.at(right_node);
            (*it).operator[](0)-= delta;
        }





        //if(debug()) std::cout<<"messages added "<<std::endl;

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
