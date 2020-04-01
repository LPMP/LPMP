#pragma once
#include <andres/graph/digraph.hxx>
#include <unordered_map>
#include <unordered_set>
#include <stack>


namespace LPMP {

template<class LDP_STRUCT>
class ldp_single_node_cut_factor
{
public:
	constexpr static std::size_t no_edge_active = std::numeric_limits<std::size_t>::infinity();

	//template<class LPD_STRUCT> ldp_single_node_cut_factor(const LPD_STRUCT& ldpStruct);
	ldp_single_node_cut_factor(const LDP_STRUCT& ldpStruct,size_t nID,bool isOut):
		baseGraph(ldpStruct.getGraph()),
		liftedGraph(ldpStruct.getGraphLifted()),
		nodeID(nID),
		ldpStructure(ldpStruct),
		isOutFlow(isOut)
	{
		numberOfEdges=baseGraph.numberOfEdgesFromVertex(nodeID);
		numberOfLiftedEdges=liftedGraph.numberOfEdgesFromVertex(nodeID);

		primal_=numberOfEdges;  //corresponds to no edge active
		optimalSolution=numberOfEdges;

		if(isOutFlow){
			minLayer=ldpStruct.getGroupIndex(nodeID);
			maxLayer=minLayer+ldpStruct.getGapLifted(); //some method that returns max time gap lifted
		}
		else{
			maxLayer=ldpStruct.getGroupIndex(nodeID);
			minLayer=std::max(0,int(maxLayer)-int(ldpStruct.getGapLifted()));
		}

		baseCosts=std::vector<double>(numberOfEdges);
		liftedCosts=std::vector<double>(numberOfLiftedEdges);
		solutionCosts=std::vector<double>(numberOfEdges+1,0);



	}

	virtual size_t getNeighborBaseEdge(size_t firstNode,size_t neighborIndex)=0;
	virtual size_t getNeighborBaseVertex(size_t firstNode,size_t neighborIndex)=0;
	virtual size_t getNeighborLiftedEdge(size_t firstNode,size_t neighborIndex)=0;
	virtual size_t getNeighborLiftedVertex(size_t firstNode,size_t neighborIndex)=0;
	virtual size_t numberOfNeighborsBase(size_t nodeIndex)=0;
	virtual size_t numberOfNeighborsLifted(size_t nodeIndex)=0;
	virtual bool inInRange(size_t nodeIndex)=0;
	virtual bool reachable(size_t firstVertex,size_t secondVertex)=0;
	virtual ~ldp_single_node_cut_factor()=0;
	virtual size_t getVertexToReach()=0;
	virtual std::pair<bool,size_t> findEdgeBase(size_t firstNode,size_t secondNode)=0;
	virtual std::pair<bool,size_t> findEdgeLifted(size_t firstNode,size_t secondNode)=0;

	//ldp_single_node_cut_factor(const std::size_t nr_outgoing_base_edges, const std::size_t nr_outgoing_lifted_edges);
	double LowerBound() const{return solutionCosts[optimalSolution]; }
	double EvaluatePrimal() const {return solutionCosts[primal_]; }


	template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
	template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

	//auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }//TODO change this. This will not work with so many variables
	auto export_variables() { return std::tie(solutionCosts); }//TODO change this. This will not work with so many variables

	//void init_primal() { primal_ = no_edge_active; }
	void init_primal() { primal_ = numberOfEdges; }

	void updateCost(const double value,const size_t index);

	const andres::graph::Digraph<>&& getBaseGraph() const {
		return baseGraph;
	}

	const LDP_STRUCT& getLdpStructure() const {
		return ldpStructure;
	}

	const andres::graph::Digraph<>&& getLiftedGraph() const {
		return liftedGraph;
	}

	std::size_t getMaxLayer() const {
		return maxLayer;
	}

	std::size_t getMinLayer() const {
		return minLayer;
	}

	std::size_t getNumberOfEdges() const {
		return numberOfEdges;
	}

	const std::size_t nodeID;


private:
	void updateValues(size_t liftedEdgeID=numberOfLiftedEdges);  //Highest number: update all.For base edge update, simpler procedure
	//TODO implement update for base edge cost


	std::size_t primal_; // the incoming resp. outgoing edge that is active.
	//std::vector<bool> primal_; // the incoming resp. outgoing edge that is active.
	std::size_t optimalSolution;
	std::size_t numberOfEdges;
	std::size_t numberOfLiftedEdges;

	// double primalValue_=0;

	std::size_t minLayer;
	std::size_t maxLayer;

	bool isOutFlow;

	const andres::graph::Digraph<>& baseGraph;
	const andres::graph::Digraph<>& liftedGraph;
	const LDP_STRUCT& ldpStructure;


	std::vector<double> baseCosts;
	std::vector<double> liftedCosts;
	std::unordered_map<size_t,double> valuesStructure;  //vertexID->value,

	std::vector<double> solutionCosts;

	size_t decodeIndex(size_t index){
		return index;  //TODO implement decoding of message index into edge/lifted edge indices
	}


};

template<class LDP_STRUCT>
class ldp_node_cut_outflow: public ldp_single_node_cut_factor<LDP_STRUCT>
{
public:
	size_t getNeighborBaseEdge(size_t firstNode,size_t neighborIndex){
		return getBaseGraph().edgeFromVertex(firstNode,neighborIndex);
	}
	size_t getNeighborBaseVertex(size_t firstNode,size_t neighborIndex){
		return getBaseGraph().vertexFromVertex(firstNode,neighborIndex);
	}
	size_t getNeighborLiftedEdge(size_t firstNode,size_t neighborIndex){
		return getLiftedGraph().edgeFromVertex(firstNode,neighborIndex);
	}
	size_t getNeighborLiftedVertex(size_t firstNode,size_t neighborIndex){
		return getLiftedGraph().vertexFromVertex(firstNode,neighborIndex);
	}
	bool inInRange(size_t nodeIndex){
		return getLdpStructure().getGroupIndex(nodeIndex)<=maxLayer;
	}
	size_t numberOfNeighborsBase(size_t nodeIndex){
		return getBaseGraph().numberOfEdgesFromVertex(nodeIndex);
	}
    size_t numberOfNeighborsLifted(size_t nodeIndex){
    	return getLiftedGraph().numberOfEdgesFromVertex(nodeIndex);
    }
    bool reachable(size_t firstVertex,size_t secondVertex){
        return getLdpStructure().isReachable(firstVertex,secondVertex);
    }
    size_t getVertexToReach(){
    	return getLdpStructure().getTerminalNode;
    }
    std::pair<bool,size_t> findEdgeBase(size_t firstNode,size_t secondNode){
    	return getBaseGraph().findEdge(firstNode,secondNode);
    }
    std::pair<bool,size_t> findEdgeLifted(size_t firstNode,size_t secondNode){
    	return getLiftedGraph().findEdge(firstNode,secondNode);
    }

	~ldp_node_cut_outflow(){ }

};


template<class LDP_STRUCT>
class ldp_node_cut_inflow: public ldp_single_node_cut_factor<LDP_STRUCT>
{
public:
	size_t getNeighborBaseEdge(size_t firstNode,size_t neighborIndex){
		return getBaseGraph().edgeToVertex(firstNode,neighborIndex);
	}
	size_t getNeighborBaseVertex(size_t firstNode,size_t neighborIndex){
		return getBaseGraph().vertexToVertex(firstNode,neighborIndex);
	}
	size_t getNeighborLiftedEdge(size_t firstNode,size_t neighborIndex){
		return getLiftedGraph().edgeToVertex(firstNode,neighborIndex);
	}
	size_t getNeighborLiftedVertex(size_t firstNode,size_t neighborIndex){
		return getLiftedGraph().vertexToVertex(firstNode,neighborIndex);
	}
	bool inInRange(size_t nodeIndex){
		return getLdpStructure().getGroupIndex(nodeIndex)>=minLayer;
	}
	size_t numberOfNeighborsBase(size_t nodeIndex){
		return getBaseGraph().numberOfEdgesToVertex(nodeIndex);
	}
	size_t numberOfNeighborsLifted(size_t nodeIndex){
		return getLiftedGraph().numberOfEdgesToVertex(nodeIndex);
	}
	bool reachable(size_t firstVertex,size_t secondVertex){
		return getLdpStructure().isReachable(secondVertex,firstVertex);
	}
	size_t getVertexToReach(){
		return getLdpStructure().getSourceNode;
	}
	std::pair<bool,size_t> findEdgeBase(size_t firstNode,size_t secondNode){
		return getBaseGraph().findEdge(secondNode,firstNode);
	}
	std::pair<bool,size_t> findEdgeLifted(size_t firstNode,size_t secondNode){
		return getLiftedGraph().findEdge(secondNode,firstNode);
	}
	~ldp_node_cut_inflow(){ }

};


//    template<class LDP_STRUCT>
//    inline ldp_single_node_cut_factor<LPD_STRUCT>::ldp_single_node_cut_factor(const LPD_STRUCT& ldpStruct,size_t nodeID):
//	baseGraph(ldpStruct.getGraph()),
//	liftedGraph(ldpStruct.getGraphLifted())
//	{
//
//
//    }

template<class LDP_STRUCT>
inline void ldp_single_node_cut_factor<LDP_STRUCT>::updateCost(const double value,const size_t index){
	size_t myIndex=decodeIndex(index);
	if(myIndex<numberOfEdges){ //update in base edge
		baseCosts[myIndex]+=value;
		solutionCosts[myIndex]+=value;
		if(myIndex==optimalSolution){
			if(value>0){
				for (int i = 0; i < numberOfEdges; ++i) {
					if(solutionCosts[i]<solutionCosts[optimalSolution]){
						optimalSolution=i;
						primal_=i;
					}
				}
			}
		}
		else if(solutionCosts[myIndex]<solutionCosts[optimalSolution]){
			optimalSolution=myIndex;
			primal_=myIndex; //maybe not

		}

	}
	else{ //update in lifted edge
		liftedCosts[myIndex-numberOfEdges]+=value;
		updateValues(myIndex-numberOfEdges);
	}
}


template<class LDP_STRUCT>
inline void ldp_single_node_cut_factor<LDP_STRUCT>::updateValues(size_t liftedEdgeID){
	std::unordered_map<size_t,size_t> indexStructure; //vertex->vertex. For reconstruction of opt. solution in lifted edges
	//In case that we do not need the reconstruction, unordered_set will be enough
	bool oneEdgeUpdate=liftedEdgeID>=numberOfLiftedEdges;
	size_t vertexToReach=vertexToReach();
	if(oneEdgeUpdate){
		//vertexToReach=liftedGraph.edgeFromVertex(nodeID,liftedEdgeID);
		vertexToReach=getNeighborLiftedEdge(nodeID,liftedEdgeID);
	}

	std::stack<size_t> nodeStack;
	nodeStack.push(nodeID);

	while(!nodeStack.empty()){
		size_t currentNode=nodeStack.top();
		bool descClosed=true;
		double minValue=0;
		size_t minValueIndex=vertexToReach();

		for (int i = 0; i < numberOfNeighborsBase(currentNode); ++i) {
			//size_t desc=baseGraph.vertexFromVertex(currentNode,i);
			size_t desc=getNeighborBaseVertex(currentNode,i);
			if(inInRange(desc)){
				if(indexStructure.count(desc)>0||(oneEdgeUpdate&&!reachable(currentNode,vertexToReach))){  //descendant closed
					if(descClosed&&valuesStructure.count(desc)>0){
						if(minValue>valuesStructure[desc]){
							minValue=valuesStructure[desc];
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
		if(descClosed){ //Close node if all descendants are closed
			if(currentNode==nodeID){

				double bestValue=0;
				size_t bestEdge=numberOfEdges;
				for (int i = 0; i < numberOfEdges; ++i) {
					size_t vertex=getNeighborBaseVertex(nodeID,i);
					double valueToAdd=0;
					if(valuesStructure.count(vertex)>0){
						valueToAdd=valuesStructure[vertex];
					}
					else{
						auto findEdge=findEdgeLifted(nodeID,vertex);
						if(findEdge.first){
							valueToAdd=liftedCosts[findEdge.second];
						}

					}

					solutionCosts[i]=baseCosts[i]+valueToAdd;
					if(solutionCosts[i]<bestValue){
						bestValue=solutionCosts;
						bestEdge=i;
					}
				}
				optimalSolution=bestEdge;
				primal_=bestEdge;//Do this or not?

			}
			else{
				auto findEdge=findEdgeLifted(nodeID,currentNode);
				if(findEdge.first){
					minValue+=liftedCosts[findEdge.second];  //add lifted edge cost if the edge exists
				}
				auto findEdgeB=findEdgeBase(nodeID,currentNode);
				if(minValue<0||(findEdgeB.first&&minValue<liftedCosts[findEdge.second])){  //store only negative values or values needed to correct solutionCosts
					valuesStructure[currentNode]=minValue;
				}
				else{
					valuesStructure.erase(currentNode);
				}
				indexStructure[currentNode]=minValueIndex; //marking the node as closed. Points to the best descendant even if not stored
				//in valuesStructure
			}
			nodeStack.pop();
		}
	}


}








class ldp_mcf_single_node_cut_message
{
	template<typename SINGLE_NODE_CUT_FACTOR>
	void RepamLeft(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
	{
		r.updateCost(msg,msg_dim);
		//Update costs in vectors, run updateValues or simpleUpdateValues

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

}
