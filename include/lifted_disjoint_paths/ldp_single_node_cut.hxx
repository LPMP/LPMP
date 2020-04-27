#pragma once
#include <andres/graph/digraph.hxx>
#include <unordered_map>
#include <unordered_set>
#include <stack>


namespace LPMP {

struct baseAndLiftedMessages{
	std::unordered_map<size_t,double> baseMessages;
	std::unordered_map<size_t,double> liftedMessages;
};

template<class LDP_INSTANCE>
class ldp_single_node_cut_factor
{
public:
	//constexpr static std::size_t no_edge_active = std::numeric_limits<std::size_t>::infinity();

	//template<class LPD_STRUCT> ldp_single_node_cut_factor(const LPD_STRUCT& ldpStruct);
	ldp_single_node_cut_factor(const LDP_INSTANCE& ldpInst,size_t nID,bool isOut):
		baseGraph(ldpInst.getGraph()),
		//liftedGraph(ldpStruct.getGraphLifted()),
		nodeID(nID),
		ldpInstance(ldpInst),
		isOutFlow(isOut),
		nodeNotActive(nID)
{


		primal_=nodeNotActive;  //corresponds to no edge active
		optimalSolution=nodeNotActive;

		if(isOutFlow){
			minLayer=ldpInst.getGroupIndex(nodeID);
			maxLayer=minLayer+ldpInst.getGapLifted(); //some method that returns max time gap lifted
		}
		else{
			maxLayer=ldpInst.getGroupIndex(nodeID);
			minLayer=std::max(0,int(maxLayer)-int(ldpInst.getGapLifted()));
		}

		baseCosts=std::unordered_map<size_t,double>();  //TODO change to maps nodeID->cost. Should be easier for messages
		liftedCosts=std::unordered_map<size_t,double>();
		solutionCosts=std::unordered_map<size_t,double>();
		solutionCosts[nodeNotActive]=0;
		vsUpToDate=false;

}


	double LowerBound() const;
	double EvaluatePrimal() const{
		double value=solutionCosts.at(primal_);
		return value;
	}
	std::unordered_map<size_t,double> adjustCostsAndSendMessages();



	//baseAndLiftedMessages adjustCostsAndSendMessagesLifted();

	//std::unordered_map<size_t,std::unordered_set<size_t>> initPredecessorsInIndexStr();


	template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
	template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

	//auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }//TODO change this. This will not work with so many variables
	auto export_variables() { return std::tie(solutionCosts); }//?What comes here?

	void init_primal() { primal_ = nodeNotActive; }

	void updateCostSimple(const double value,const size_t vertexIndex,bool isLifted);

	void updateValues() const;

	const andres::graph::Digraph<>& getBaseGraph() const {
		return baseGraph;
	}

	const std::size_t nodeID;
	const size_t nodeNotActive;


private:

	size_t getNeighborBaseVertex(size_t firstNode,size_t neighborIndex) const{
		if(isOutFlow){
			return baseGraph.vertexFromVertex(firstNode,neighborIndex);
		}
		else{
			return baseGraph.vertexToVertex(firstNode,neighborIndex);
		}
	}
	size_t numberOfNeighborsBase(size_t nodeIndex) const {
		if(isOutFlow){
			return baseGraph.numberOfEdgesFromVertex(nodeIndex);
		}
		else{
			return baseGraph.numberOfEdgesToVertex(nodeIndex);
		}
	}
	bool isInRange(size_t nodeIndex) const {
		if(isOutFlow){
			return ldpInstance.getGroupIndex(nodeIndex)<=maxLayer;
		}
		else{
			return ldpInstance.getGroupIndex(nodeIndex)>=minLayer;
		}
	}



	std::size_t primal_; // the incoming resp. outgoing edge that is active.
	mutable std::size_t optimalSolution;

	std::size_t minLayer;
	std::size_t maxLayer;

	const bool isOutFlow;

	const andres::graph::Digraph<>& baseGraph;
	//const andres::graph::Digraph<>& liftedGraph;
	const LDP_INSTANCE& ldpInstance;


	mutable std::unordered_map<size_t,double> baseCosts;
	std::unordered_map<size_t,double> liftedCosts;
	mutable std::unordered_map<size_t,double> solutionCosts;
	mutable std::unordered_map<size_t,std::unordered_set<size_t>> indexStructure;

	mutable std::unordered_map<size_t,double> valuesStructure;  //For DFS procedure

	mutable bool vsUpToDate;

	//	 size_t getNeighborBaseEdge(size_t firstNode,size_t neighborIndex){
	//		 if(isOutFlow){
	//			 return baseGraph.edgeFromVertex(firstNode,neighborIndex);
	//		 }
	//		 else{
	//			 return baseGraph.edgeToVertex(firstNode,neighborIndex);
	//		 }
	//	 }
	//	 size_t getNeighborLiftedEdge(size_t firstNode,size_t neighborIndex){
	//		 if(isOutFlow){
	//			 return liftedGraph.edgeFromVertex(firstNode,neighborIndex);
	//		 }
	//		 else{
	//			 return liftedGraph.edgeToVertex(firstNode,neighborIndex);
	//		 }
	//	 }
	//	 size_t getNeighborLiftedVertex(size_t firstNode,size_t neighborIndex){
	//		 if(isOutFlow){
	//			 return liftedGraph.vertexFromVertex(firstNode,neighborIndex);
	//		 }
	//		 else{
	//			 return liftedGraph.vertexToVertex(firstNode,neighborIndex);
	//		 }
	//	 }

	//	 size_t numberOfNeighborsLifted(size_t nodeIndex){
	//		 if(isOutFlow){
	//			 return liftedGraph.numberOfEdgesFromVertex(nodeIndex);
	//		 }
	//		 else{
	//			 return liftedGraph.numberOfEdgesToVertex(nodeIndex);
	//		 }
	//	 }

	//	 bool reachable(size_t firstVertex,size_t secondVertex){
	//		 if(isOutFlow){
	//			 return ldpStructure.isReachable(firstVertex,secondVertex);
	//		 }
	//		 else{
	//			 return ldpStructure.isReachable(secondVertex,firstVertex);
	//		 }
	//	 }


	//	 std::pair<bool,size_t> findEdgeBase(size_t firstNode,size_t secondNode){
	//		 if(isOutFlow){
	//			 return baseGraph.findEdge(firstNode,secondNode);
	//		 }
	//		 else{
	//			 return baseGraph.findEdge(secondNode,firstNode);
	//		 }
	//
	//	 }
	//	 std::pair<bool,size_t> findEdgeLifted(size_t firstNode,size_t secondNode){
	//		 if(isOutFlow){
	//			 return liftedGraph.findEdge(firstNode,secondNode);
	//		 }
	//		 else{
	//			 return liftedGraph.findEdge(secondNode,firstNode);
	//		 }
	//	 }


		size_t getVertexToReach(){
			if(isOutFlow){
				return ldpInstance.getTerminalNode;
			}
			else{
				return ldpInstance.getSourceNode;
			}
		}

};



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


template<class LDP_INSTANCE>
inline std::unordered_map<size_t,double> ldp_single_node_cut_factor<LDP_INSTANCE>::adjustCostsAndSendMessages(){
	if(!vsUpToDate) updateValues();
	double minValue=solutionCosts[optimalSolution];
	std::unordered_map<size_t,double> messages;
	for (auto it=solutionCosts.begin();it!=solutionCosts.end();it++) {
		if(it->first==optimalSolution) continue;
		double delta=it->second-minValue;
		messages[it->first]=delta;
		baseCosts[it->first]-=delta;
		it->second=minValue;
	}
	return messages;
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateCostSimple(const double value,const size_t vertexIndex,bool isLifted){//Only cost change
	if(!isLifted){ //update in base edge
		baseCosts[vertexIndex]+=value;
		solutionCosts[vertexIndex]+=value;
	}
	else{ //update in lifted edge
		liftedCosts[vertexIndex]+=value;
		valuesStructure[vertexIndex]+=value;
		vsUpToDate=false;
	}
}


template<class LDP_INSTANCE>
inline double ldp_single_node_cut_factor<LDP_INSTANCE>::LowerBound() const{//TODO store info about how valuesStructures changed. At least max time layer of changed lifted edge
	if(!vsUpToDate){
		updateValues();
		return solutionCosts.at(optimalSolution);
	}
	else{
		double minValue=solutionCosts.at(optimalSolution);
		size_t minVertex=optimalSolution;
		for (std::pair<size_t,double> it :solutionCosts) {
			if(it.second<minValue){
				minValue=it.second;
				minVertex=it.first;
			}
		}
		optimalSolution=minVertex;
		return solutionCosts.at(optimalSolution);
	}
}


template<class LDP_INSTANCE>
inline void ldp_single_node_cut_factor<LDP_INSTANCE>::updateValues() const{
	//TODO check if everything works fine for edge from v to t
	std::unordered_set<size_t> closedVertices;
	indexStructure=std::unordered_map<size_t,std::unordered_set<size_t>>();
	valuesStructure=std::unordered_map<size_t,double>();
	//size_t vertexToReach=getVertexToReach();

	std::stack<size_t> nodeStack;
	nodeStack.push(nodeID);

	while(!nodeStack.empty()){
		size_t currentNode=nodeStack.top();
		bool descClosed=true;
		double minValue=0;
		std::unordered_set<size_t> minValueIndices;
		//size_t minValueIndex=vertexToReach;

		for (int i = 0; i < numberOfNeighborsBase(currentNode); ++i) {
			//size_t desc=baseGraph.vertexFromVertex(currentNode,i);
			size_t desc=getNeighborBaseVertex(currentNode,i);
			if(isInRange(desc)){
				//if(indexStructure.count(desc)>0||(oneEdgeUpdate&&!reachable(currentNode,vertexToReach))){  //descendant closed
				if(closedVertices.count(desc)>0){  //descendant closed
					if(descClosed&&valuesStructure.count(desc)>0){
						if(minValue>=valuesStructure[desc]&&valuesStructure[desc]<0){
							minValue=valuesStructure[desc];
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
			if(currentNode==nodeID){  //all nodes closed, compute solution values
				double bestValue=0;
				size_t bestVertex=nodeNotActive;
				for (auto it=baseCosts.begin();it!=baseCosts.end();it++) {
					size_t vertex=it->first;
					double baseCost=it->second;
					double valueToAdd=0;
					if(valuesStructure.count(vertex)>0){
						valueToAdd=valuesStructure[vertex];
					}
					else{
						if(liftedCosts.count(vertex)>0){
							valueToAdd=liftedCosts.at(vertex);
						}
					}
					solutionCosts[vertex]=baseCost+valueToAdd;
					if(solutionCosts[vertex]<bestValue){
						bestValue=solutionCosts[vertex];
						bestVertex=vertex;
					}
				}
				optimalSolution=bestVertex;
//				for (int i = 0; i < numberOfNeighborsBase(nodeID); ++i) {//Only for index structure
//					size_t desc=getNeighborBaseVertex(nodeID,i);
//					if(isInRange(desc)){
//						if(solutionCosts[desc]==bestValue){
//							indexStructure[nodeID].insert(desc);
//						}
//					}
//				}

			}
			else{
				//if(minValue<0){//This for cycle is only for index structure creation
//				for (int i = 0; i < numberOfNeighborsBase(currentNode); ++i) {
//					size_t desc=getNeighborBaseVertex(currentNode,i);
//					if(isInRange(desc)){
//						if(isInRange(desc)){
//							if(valuesStructure.count(desc)>0&&valuesStructure[desc]==minValue){
//								indexStructure[currentNode].insert(desc);
//							}
//						}
//					}
//				}
				//}
				if(liftedCosts.count(currentNode)>0){
					minValue+=liftedCosts.at(currentNode);  //add lifted edge cost if the edge exists
				}
				if(minValue<0||(baseCosts.count(currentNode)>0&&minValue<liftedCosts.at(currentNode))){  //store only negative values or values needed to correct solutionCosts
					valuesStructure[currentNode]=minValue;
				}
				else{
					valuesStructure.erase(currentNode);
				}
				closedVertices.insert(currentNode); //marking the node as closed.
			}
			nodeStack.pop();
		}
	}

	vsUpToDate=true;

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

}
