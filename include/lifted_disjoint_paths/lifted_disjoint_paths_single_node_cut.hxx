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
		vsUpToDate=false;



	}


	//ldp_single_node_cut_factor(const std::size_t nr_outgoing_base_edges, const std::size_t nr_outgoing_lifted_edges);
	double LowerBound() {//TODO store info about how valuesStructures changed. At least max time layer of changed lifted edge
		if(!vsUpToDate){
			updateValues();
		}
		else{
			double minValue=solutionCosts[optimalSolution];
			size_t minIndex=optimalSolution;
			for (int i = 0; i < numberOfEdges; ++i) {
				if(solutionCosts[i]<minValue){
					minValue=solutionCosts[i];
					minIndex=i;
				}
			}
			optimalSolution=minIndex;

		}
		return solutionCosts[optimalSolution];
	}
	double EvaluatePrimal() const {return solutionCosts[primal_]; }


	template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
	template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

	//auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }//TODO change this. This will not work with so many variables
	auto export_variables() { return std::tie(solutionCosts); }//TODO change this. This will not work with so many variables

	//void init_primal() { primal_ = no_edge_active; }
	void init_primal() { primal_ = numberOfEdges; }

	void updateCostFull(const double value,const size_t index);
	void updateCostSimple(const double value,const size_t index);

	const andres::graph::Digraph<>& getBaseGraph() const {
		return baseGraph;
	}



	const std::size_t nodeID;
	virtual ~ldp_single_node_cut_factor()=0;


private:
	 size_t getNeighborBaseEdge(size_t firstNode,size_t neighborIndex){
		 if(isOutFlow){
			 return baseGraph.edgeFromVertex(firstNode,neighborIndex);
		 }
		 else{
			 return baseGraph.edgeToVertex(firstNode,neighborIndex);
		 }
	 }
	 size_t getNeighborBaseVertex(size_t firstNode,size_t neighborIndex){
		 if(isOutFlow){
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
	 size_t numberOfNeighborsBase(size_t nodeIndex){
		 if(isOutFlow){
			 return baseGraph.numberOfEdgesFromVertex(nodeIndex);
		 }
		 else{
			 return baseGraph.numberOfEdgesToVertex(nodeIndex);
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
	 bool inInRange(size_t nodeIndex){
		 if(isOutFlow){
			 return ldpStructure.getGroupIndex(nodeIndex)<=maxLayer;
		 }
		 else{
			 return ldpStructure.getGroupIndex(nodeIndex)>=minLayer;
		 }
	 }
	 bool reachable(size_t firstVertex,size_t secondVertex){
		 if(isOutFlow){
			 return ldpStructure.isReachable(firstVertex,secondVertex);
		 }
		 else{
			 return ldpStructure.isReachable(secondVertex,firstVertex);
		 }
	 }

	 size_t getVertexToReach(){
		 if(isOutFlow){
			 return ldpStructure.getTerminalNode;
		 }
		 else{
			 return ldpStructure.getSourceNode;
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
	 std::pair<bool,size_t> findEdgeLifted(size_t firstNode,size_t secondNode){
		 if(isOutFlow){
			 return liftedGraph.findEdge(firstNode,secondNode);
		 }
		 else{
			 return liftedGraph.findEdge(secondNode,firstNode);
		 }
	 }




	void updateValues(size_t liftedEdgeID=std::numeric_limits<std::size_t>::infinity());  //Highest number: update all.For base edge update, simpler procedure
	//TODO implement update for base edge cost


	std::size_t primal_; // the incoming resp. outgoing edge that is active.
	//std::vector<bool> primal_; // the incoming resp. outgoing edge that is active.
	std::size_t optimalSolution;
	std::size_t numberOfEdges;
	std::size_t numberOfLiftedEdges;

	// double primalValue_=0;

	std::size_t minLayer;
	std::size_t maxLayer;

	const bool isOutFlow;

	const andres::graph::Digraph<>& baseGraph;
	const andres::graph::Digraph<>& liftedGraph;
	const LDP_STRUCT& ldpStructure;


	std::vector<double> baseCosts;
	std::vector<double> liftedCosts;
	std::unordered_map<size_t,double> valuesStructure;  //vertexID->value,

	std::vector<double> solutionCosts;

	bool vsUpToDate;

	size_t decodeIndex(size_t index){
		return index;  //TODO implement decoding of message index into edge/lifted edge indices
	}


};


template<class LDP_STRUCT>
inline void ldp_single_node_cut_factor<LDP_STRUCT>::updateCostFull(const double value,const size_t index){//Includes DFS update
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
inline void ldp_single_node_cut_factor<LDP_STRUCT>::updateCostSimple(const double value,const size_t index){//Only cost change
	size_t myIndex=decodeIndex(index);  //global variable index to factor variable index
	if(myIndex<numberOfEdges){ //update in base edge
		baseCosts[myIndex]+=value;
		solutionCosts[myIndex]+=value;
	}
	else{ //update in lifted edge
		liftedCosts[myIndex-numberOfEdges]+=value;
		size_t vertex=getNeighborLiftedVertex(nodeID,myIndex-numberOfEdges);
		valuesStructure[vertex]+=value;
		vsUpToDate=false;
	}
}



template<class LDP_STRUCT>
inline void ldp_single_node_cut_factor<LDP_STRUCT>::updateValues(size_t liftedEdgeID){
	//TODO probably switch from liftedEdgeID to the index of the max layer that has been changed
	std::unordered_map<size_t,size_t> indexStructure; //vertex->vertex. For reconstruction of opt. solution in lifted edges
	//In case that we do not need the reconstruction, unordered_set will be enough
	bool oneEdgeUpdate=liftedEdgeID>=numberOfLiftedEdges;
	size_t vertexToReach=getVertexToReach();
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
		size_t minValueIndex=vertexToReach;

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

	vsUpToDate=true;

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

}
