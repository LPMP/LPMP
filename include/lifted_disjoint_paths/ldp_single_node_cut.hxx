#pragma once
#include <andres/graph/digraph.hxx>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <array>
#include <list>
#include <set>
#include <config.hxx>
#include<lifted_disjoint_paths/ldp_functions.hxx>

namespace LPMP {




struct StrForTopDownUpdate{

    StrForTopDownUpdate(const std::vector<double>& bCosts,const std::vector<double>& lCosts,const size_t centralNodeID,const size_t mostDistantNodeID,const size_t vertexToReach):
	baseCosts(bCosts),
    liftedCosts(lCosts),
    solutionCosts(bCosts.size()+1),
    //solutionCosts(solCosts),
	optValue(0),
    nodeID(centralNodeID),
    topDownValuesStructure(centralNodeID,mostDistantNodeID,0),
    topDownVertexIDStructure(centralNodeID,mostDistantNodeID,vertexToReach),
     optBaseIndex(0) //maybe set to node not active
	{

	}

    const size_t nodeID;
    const std::vector<double>& baseCosts;
    const std::vector<double>& liftedCosts;

    std::vector<double> solutionCosts;
    ShiftedVector<double> topDownValuesStructure;
    ShiftedVector<size_t> topDownVertexIDStructure;

    size_t optBaseIndex;
    double optValue;




};

template<class LDP_INSTANCE>
class ldp_single_node_cut_factor
{
public:
    ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut);

	void initBaseCosts(double fractionBase);

	void initLiftedCosts(double fractionLifted);

	const std::vector<double>& getLiftedCosts() const {
		return liftedCosts;
    }

	const std::vector<double>& getBaseCosts() const {
		return baseCosts;
	}

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


    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(baseCosts, liftedCosts); }

    auto export_variables() { return std::tie(baseCosts, liftedCosts); }

	void updateCostSimple(const double value,const size_t vertexIndex,bool isLifted);
	void updateNodeCost(const double value);


    double getOneBaseEdgeMinMarginal(size_t vertex) const;
    std::vector<double> getAllBaseMinMarginals()const;
    std::vector<double> getAllBaseMinMarginalsForMCF() const;
    std::vector<double> getAllLiftedMinMarginals(std::vector<double>* pLocalBaseCosts=nullptr) const;
	double getNodeMinMarginal()const;
    double getOneLiftedMinMarginal(size_t indexOfLiftedEdge)const;


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
  //  void updateValues() const;
    //methods used in computing lower bound or min marginals
    ShiftedVector<bool> topDownUpdate(StrForTopDownUpdate& myStr, const size_t vertexIDToIgnore) const;
    ShiftedVector<bool> topDownUpdate(StrForTopDownUpdate& myStr) const;
    std::unordered_map<size_t,double> bottomUpUpdate(const StrForTopDownUpdate& myStr, const size_t vertex, ShiftedVector<size_t>& indexStr, ShiftedVector<double> &pBUValuesStr, const ShiftedVector<bool> &verticesInScope, ShiftedVector<bool>* pClosedVert=nullptr)const;
    void updateOptimal() const;

    //Obtain IDs of vertices of lifted edges that are part of a found optimal solution
    std::list<size_t> getOptLiftedFromIndexStr(const StrForTopDownUpdate& myStr)const;


    //Methods for exploring graph structures
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



    size_t primalBase_;
    std::unordered_set<size_t> primalLifted_;

    mutable size_t optimalSolutionBase; //Contains vertices w.r.t. their order for central vertex

    std::size_t mostDistantNeighborID;

	const bool isOutFlow;

	const andres::graph::Digraph<>& baseGraph;
	const andres::graph::Digraph<>& liftedGraph;
	const LDP_INSTANCE& ldpInstance;

    std::vector<double> baseCosts;
	std::vector<double> liftedCosts;
    double nodeCost;
	mutable std::vector<double> solutionCosts;



    std::unordered_map<size_t,size_t> baseIDToIndex;
	std::unordered_map<size_t,size_t> liftedIDToOrder;


	mutable bool optLiftedUpToDate;
	mutable bool optBaseUpToDate;

    mutable bool solutionCostsUpToDate;
    mutable bool optValueUpToDate;

	mutable double optValue;


	std::vector<size_t> baseIDs;
	std::vector<size_t> liftedIDs;

	size_t nodeNotActive;


};


template<class LDP_INSTANCE>
inline  ldp_single_node_cut_factor<LDP_INSTANCE>::ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut):
baseGraph(ldpInst.getGraph()),
liftedGraph(ldpInst.getGraphLifted()),
nodeID(nID),
ldpInstance(ldpInst),
isOutFlow(isOut)

{


    mostDistantNeighborID=nodeID;
    initBaseCosts(0);
    initLiftedCosts(0);

	nodeNotActive=baseCosts.size();

	optLiftedUpToDate=false;
	optBaseUpToDate=false;

	primalBase_=nodeNotActive;  //corresponds to no edge active
	optimalSolutionBase=nodeNotActive;
	nodeCost=0;
	optValue=0;

    optValueUpToDate=false;
    solutionCostsUpToDate=false;



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
                double toAdd=myStr.liftedCosts.at(liftedIDToOrder.at(vertexInOptimalPath));
                optValueComputed+=toAdd;

            }
            vertexInOptimalPath=myStr.topDownVertexIDStructure[vertexInOptimalPath];
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
            for(size_t i=0;i<solutionCosts.size();i++){
                if(optValue>solutionCosts[i]){
                    optValue=solutionCosts[i];
                    optimalSolutionBase=i;
                }
            }
        }
        else{
            StrForTopDownUpdate myStr(baseCosts,liftedCosts,nodeID,mostDistantNeighborID,getVertexToReach());
            topDownUpdate(myStr);
            optValue=myStr.optValue;
            solutionCosts=myStr.solutionCosts;
            optimalSolutionBase=myStr.optBaseIndex;
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
    baseIDToIndex.clear();
	if(fractionBase==0){
		for (int i = 0; i < numberOfNeighborsBase(nodeID); ++i) {
			size_t neighborID=getNeighborBaseVertex(nodeID,i);
            if(neighborID!=getVertexToReach()){
                mostDistantNeighborID=getMoreDistantNode(mostDistantNeighborID,neighborID);
            }
			baseCosts.push_back(0);
			baseIDs.push_back(neighborID);
            baseIDToIndex[neighborID]=i;
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
            baseCosts.push_back(fractionBase*cost);
			baseIDs.push_back(neighborID);
            baseIDToIndex[neighborID]=i;

		}

	}
    solutionCosts=std::vector<double>(baseCosts.size()+1);
    optBaseUpToDate=false;
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
inline ShiftedVector<bool> ldp_single_node_cut_factor<LDP_INSTANCE>::topDownUpdate(StrForTopDownUpdate& myStr) const{
    return topDownUpdate(myStr,getVertexToReach());
}


template<class LDP_INSTANCE>
inline ShiftedVector<bool> ldp_single_node_cut_factor<LDP_INSTANCE>::topDownUpdate(StrForTopDownUpdate& myStr,const size_t vertexIDToIgnore) const{

    ShiftedVector<bool> closedVertices(nodeID,mostDistantNeighborID);


    bool vertexToIgnoreSet=false;
    size_t lastVertex=mostDistantNeighborID;
    if(vertexIDToIgnore!=getVertexToReach()){
        vertexToIgnoreSet=true;
        lastVertex=vertexIDToIgnore;
	}


    if(vertexToIgnoreSet){
        myStr.topDownValuesStructure.fillWith(0,nodeID,vertexIDToIgnore);
    }

    for (int i = 0; i < liftedIDs.size(); ++i) {
        if(!vertexToIgnoreSet||isInGivenInterval(liftedIDs.at(i),lastVertex)){
            myStr.topDownValuesStructure[liftedIDs.at(i)]=myStr.liftedCosts.at(i);
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
       //     std::cout<<"current node "<<currentNode<<std::endl;

			bool descClosed=true;
			double minValue=0;
			size_t minValueVertexID=getVertexToReach();

			for (int i = 0; i < numberOfNeighborsBase(currentNode); ++i) {
				size_t desc=getNeighborBaseVertex(currentNode,i);
                if(desc==vertexIDToIgnore||desc==getVertexToReach()) continue;

                if(isInGivenInterval(desc,mostDistantNeighborID)){

                    assert(desc<baseGraph.numberOfVertices());
                    if(closedVertices.getValue(desc)||(vertexToIgnoreSet&&!isInGivenInterval(desc, vertexIDToIgnore))){  //descendant closed
                        if(descClosed){
                            double value=myStr.topDownValuesStructure[desc];
                            if(minValue>value){
                                minValue=value;
                                minValueVertexID=desc;
                            }
                        }
                    }
					else {  //descendant not closed
                       // std::cout<<"push desc "<<desc<<std::endl;
						nodeStack.push(desc);
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
                        double baseCost=myStr.baseCosts[i];
                        if(!vertexToIgnoreSet||baseIDs[i]!=vertexIDToIgnore){

                            double valueToAdd=0;
                            if(baseIDs[i]!=getVertexToReach()){
                                valueToAdd=myStr.topDownValuesStructure[baseIDs[i]];
                            }
                            double value=baseCost+nodeCost+valueToAdd;

                            myStr.solutionCosts.at(i)=value;
                            if(value<bestValue){
                                bestValue=value;
                                bestIndex=i;

                                //std::cout<<"best vertex"<<std::endl;
                            }
                        }
						//std::cout<<"end for"<<std::endl;liftedCosts
					}

                    myStr.optBaseIndex=bestIndex;
                   // std::cout<<"best index "<<bestIndex;
//                    if(bestIndex!=nodeNotActive){
//                        std::cout<<"~"<<baseIDs.at(bestIndex)<<std::endl;
//                    }

					myStr.optValue=bestValue;
                 //   std::cout<<"best value "<<bestValue<<std::endl;
				}
				else{

					double valueToStore=minValue;

                   // std::cout<<"best desc "<<minValueVertexID<<std::endl;
                    assert(!vertexToIgnoreSet||minValueVertexID!=vertexIDToIgnore);
                    assert(currentNode<baseGraph.numberOfVertices());
                    myStr.topDownVertexIDStructure[currentNode]=minValueVertexID;
                    myStr.topDownValuesStructure[currentNode]+=minValue;

                    closedVertices.setValue(currentNode,true); //marking the node as closed.

				}
				nodeStack.pop();

			}
		}

	}
 //   std::cout<<"update values finished"<<std::endl;

    return closedVertices;

}



template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::getOneLiftedMinMarginal(size_t indexOfLiftedEdge)const{
    assert(indexOfLiftedEdge<liftedCosts.size());


    StrForTopDownUpdate strForUpdateValues(baseCosts,liftedCosts,nodeID,mostDistantNeighborID,getVertexToReach());
    ShiftedVector<bool> verticesInScope= topDownUpdate(strForUpdateValues);
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

        topDownUpdate(strForUpdateValues,liftedIDs.at(indexOfLiftedEdge));
        double restrictedOptValue=strForUpdateValues.optValue;

        double valueToReturn=origOptValue-restrictedOptValue;
		if(debug()&&valueToReturn>1e-8){
			std::cout<<"wrong min marginal"<<nodeID<<", "<<indexOfLiftedEdge<<", opt vertex but positive value "<<valueToReturn<<std::endl;
		}
		return valueToReturn;



	}
	else{

        ShiftedVector<size_t> indexStr(nodeID,mostDistantNeighborID,getVertexToReach());
        ShiftedVector<double> buValuesStructure(nodeID,mostDistantNeighborID,0);
        for (int i = 0; i < baseIDs.size(); ++i) {
            if(baseIDs.at(i)==getVertexToReach()) continue;
            indexStr[baseIDs.at(i)]=nodeID;
            buValuesStructure[baseIDs.at(i)]=baseCosts.at(i);
        }

        std::unordered_map<size_t,double> message=bottomUpUpdate(strForUpdateValues,liftedIDs[indexOfLiftedEdge],indexStr,buValuesStructure,verticesInScope);

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
inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::bottomUpUpdate(const StrForTopDownUpdate& myStr,const size_t vertex,ShiftedVector<size_t>& indexStr,ShiftedVector<double>& buValuesStr,const ShiftedVector<bool>& verticesInScope, ShiftedVector<bool>* pClosedVert)const{
  //  std::cout<<"bottom up update from "<<vertex<<std::endl;
    bool onlyOne=pClosedVert==nullptr;


	std::unordered_map<size_t,double> messages;
	if(onlyOne){
        pClosedVert=new ShiftedVector<bool> (nodeID,mostDistantNeighborID);
    }

    ShiftedVector<bool>& closedVertices=*pClosedVert;

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
                if(pred==nodeID||!verticesInScope.isWithinBounds(pred)||!verticesInScope.getValue(pred)) continue;
                if(!closedVertices.getValue(pred)){
					predClosed=false;
					myStack.push(pred);
				}
			}
            if(predClosed){

                size_t bestIndex=getVertexToReach();
                double bestValue=std::numeric_limits<double>::infinity();
                if(indexStr[currentVertex]==nodeID){
                    bestValue=buValuesStr[currentVertex]+nodeCost;
                    bestIndex=nodeID;
                }
				for (int i = 0; i < numberOfNeighborsBaseRev(currentVertex); ++i) {
                    size_t pred=getNeighborBaseVertexRev(currentVertex,i);
                    if(pred==nodeID||!verticesInScope.isWithinBounds(pred)||!verticesInScope.getValue(pred)) continue;
                    assert(closedVertices.getValue(pred));
                    double value=buValuesStr[pred];
                    if(value<bestValue){
                        bestValue=value;
                        bestIndex=pred;  //TODO check
                    }

                }

				indexStr[currentVertex]=bestIndex;
                if(isLiftedVertex(currentVertex)){
					if(onlyOne&&currentVertex!=vertex){
                        bestValue+=myStr.liftedCosts.at(liftedIDToOrder.at(currentVertex));
					}
					else{
//
                        double topDownValueOfDesc=0;

                        size_t bestDesc=myStr.topDownVertexIDStructure[currentVertex];
                        if(bestDesc!=getVertexToReach()){
                            topDownValueOfDesc=myStr.topDownValuesStructure[bestDesc];
                        }
                        double topDownCurrent=myStr.topDownValuesStructure[currentVertex];

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
        pClosedVert=nullptr;

	}
	return messages;


}


template<class LDP_INSTANCE>
inline std::vector<double> ldp_single_node_cut_factor<LDP_INSTANCE>::getAllLiftedMinMarginals(std::vector<double>* pLocalBaseCosts) const{

	std::unordered_map<size_t,double> liftedMessages;


	std::vector<double> localLiftedCosts=liftedCosts;
    std::vector<double> localBaseCosts;
    if(pLocalBaseCosts==nullptr){
        localBaseCosts=baseCosts;
    }
    else{
        localBaseCosts=*pLocalBaseCosts;
    }

    StrForTopDownUpdate myStr(localBaseCosts,localLiftedCosts,nodeID,mostDistantNeighborID,getVertexToReach());
    topDownUpdate(myStr);
   // std::cout<<"all lifted min marginals init "<<nodeID<<std::endl;

    std::list<size_t> isNotZeroInOpt=getOptLiftedFromIndexStr(myStr);
    std::unordered_set<size_t> isOneInOpt(isNotZeroInOpt.begin(),isNotZeroInOpt.end());


    double origOptValue=myStr.optValue;
    double minMarginalsImproving=0;
    double currentOptValue=myStr.optValue;


    auto listIt=isNotZeroInOpt.begin();
	while(!isNotZeroInOpt.empty()){
		size_t vertexToClose=*listIt;
//        if(debug()) std::cout<<"process vertex "<<vertexToClose<<std::endl;

        topDownUpdate(myStr,vertexToClose);
		double newOpt=myStr.optValue;
       // std::cout<<"second best, ignoring "<<vertexToClose<<std::endl;

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
                liftedMessages[*sbIt]=0;
				sbIt++;
			}
			else if(reachable(*listIt,*sbIt)){

				listIt=isNotZeroInOpt.erase(listIt);
			}
			else{

				listIt=isNotZeroInOpt.erase(listIt);
				isOneInOpt.insert(*sbIt);
                liftedMessages[*sbIt]=0;
				sbIt++;
			}


		}



		isNotZeroInOpt.erase(listIt,isNotZeroInOpt.end());
		while(sbIt!=secondBest.end()){
			isOneInOpt.insert(*sbIt);
            liftedMessages[*sbIt]=0;
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

   // std::cout<<"min marginals improve "<<minMarginalsImproving<<std::endl;
    myStr.topDownValuesStructure.fillWith(0);  //TODO implement this in reset in StrForUpdate
    myStr.topDownVertexIDStructure.fillWith(getVertexToReach());
    ShiftedVector<bool> verticesInScope= topDownUpdate(myStr);

    assert(currentOptValue+minMarginalsImproving-origOptValue>-eps);

   // std::cout<<"computed opt "<<myStr.optValue<<std::endl;
    assert(std::abs(currentOptValue-myStr.optValue)<eps);

    ShiftedVector<double> buValuesStructure(nodeID,mostDistantNeighborID,std::numeric_limits<double>::max());
    ShiftedVector<bool> closedVertices(nodeID,mostDistantNeighborID,false);

    for(size_t optVertex:isOneInOpt){
        assert(optVertex<baseGraph.numberOfVertices());
        size_t bestDesc=myStr.topDownVertexIDStructure[optVertex];
        double toSubtract=0;
        if(bestDesc!=getVertexToReach()) toSubtract=myStr.topDownValuesStructure[bestDesc];
        buValuesStructure[optVertex]=currentOptValue-toSubtract;
        closedVertices.setValue(optVertex,true);
	}


    ShiftedVector<size_t> indexStr(nodeID,mostDistantNeighborID,getVertexToReach());  //Note that vertices closed in previous for do not have valid indexStr entries
    for (int i = 0; i < baseIDs.size(); ++i) {
        if(baseIDs.at(i)==getVertexToReach()) continue;
        if(closedVertices.getValue(baseIDs.at(i))) continue;
        indexStr[baseIDs.at(i)]=nodeID;
        buValuesStructure[baseIDs.at(i)]=baseCosts.at(i);
    }

	for(size_t vertexID:liftedIDs){
        if(!closedVertices.getValue(vertexID)){
            std::unordered_map<size_t,double> newMessages=bottomUpUpdate(myStr,vertexID,indexStr,buValuesStructure,verticesInScope,&closedVertices);
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
        StrForTopDownUpdate myStr2(localBaseCosts,localLiftedCosts,nodeID,mostDistantNeighborID,getVertexToReach());
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
