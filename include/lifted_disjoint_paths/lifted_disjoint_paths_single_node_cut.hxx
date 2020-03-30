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

            //ldp_single_node_cut_factor(const std::size_t nr_outgoing_base_edges, const std::size_t nr_outgoing_lifted_edges);
            double LowerBound() const{return solutionCosts[optimalSolution]; }
            double EvaluatePrimal() const {return solutionCosts[primal_]; }


            template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(primal_); }
            template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

            auto export_variables() { return std::tie(*static_cast<std::size_t>(this)); }//TODO change this. This will not work with so many variables

            //void init_primal() { primal_ = no_edge_active; }
            void init_primal() { primal_ = numberOfEdges; }

        private:
            void updateValues(size_t liftedEdgeID=numberOfLiftedEdges,double update=0);  //Highest number: update all.For base edge update, simpler procedure

            std::size_t nodeID;
            std::size_t primal_; // the incoming resp. outgoing edge that is active.
            //std::vector<bool> primal_; // the incoming resp. outgoing edge that is active.
            std::size_t optimalSolution;
            std::size_t numberOfEdges;
            std::size_t numberOfLiftedEdges;

            double primalValue_=0;

            std::size_t minLayer;
            std::size_t maxLayer;

            bool isOutFlow;

            const andres::graph::Digraph<>& baseGraph;
            const andres::graph::Digraph<>& liftedGraph;
            const LDP_STRUCT& ldpStructure;

            std::vector<double> baseCosts;
            std::vector<double> liftedCosts;
            std::unordered_map<size_t,double> valuesStructure;  //vertexID->value,
            std::unordered_map<size_t,size_t> indexStructure; //vertex->vertex. For reconstruction of opt. solution in lifted edges

            std::vector<double> solutionCosts;


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
        inline void ldp_single_node_cut_factor<LDP_STRUCT>::updateValues(size_t liftedEdgeID,double update){
        bool oneEdge=liftedEdgeID>=numberOfLiftedEdges;
        size_t vertexToReach=ldpStructure.getTerminalNode();
        if(oneEdge){
        	vertexToReach=liftedGraph.edgeFromVertex(nodeID,liftedEdgeID);
        }

    	std::stack<size_t> nodeStack;
        	//std::stack<size_t> parentStack;
        	nodeStack.push(nodeID);

        	while(!nodeStack.empty()){
        		size_t currentNode=nodeStack.top();
        		bool descClosed=true;
        		double minValue=0;
        		size_t minValueIndex=ldpStructure.getTerminalNode;

        		for (int i = 0; i < baseGraph.numberOfEdgesFromVertex(currentNode); ++i) {
					size_t desc=baseGraph.vertexFromVertex(currentNode,i);
					if(ldpStructure.getGroupIndex(desc)<=maxLayer&&(!oneEdge||ldpStructure.isReachable(currentNode,vertexToReach))){
						if(indexStructure.count(desc)>0){  //descendant closed
							if(descClosed&&valuesStructure.count(desc)>0){
								if(minValue>valuesStructure[desc]){
									minValue=valuesStructure[desc];
									minValueIndex=desc;
								}
							}
						}
						else{  //descendant not closed
							nodeStack.push(desc);
							//parentStack.push(currentNode);
							descClosed=false;
						}
					}
				}
        		if(descClosed){ //Close node if all descendants are closed
        			if(currentNode==nodeID){

        				double bestValue=0;
        				size_t bestEdge=numberOfEdges;
        				for (int i = 0; i < numberOfEdges; ++i) {
        					size_t vertex=baseGraph.vertexFromVertex(nodeID,i);
        					double valueToAdd=0;
        					if(valuesStructure.count(vertex)>0){
        						valueToAdd=valuesStructure[vertex];
        					}
        					else{
        						auto findEdge=liftedGraph.findEdge(nodeID,vertex);
        						if(findEdge.first){
        							valueToAdd=liftedCosts[findEdge.second];
        						}
        						size_t forwardVertex=indexStructure[vertex];
        						if(forwardVertex!=ldpStructure.getTerminalNode){
        							valueToAdd+=valuesStructure[forwardVertex];
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
        				auto findEdge=liftedGraph.findEdge(nodeID,currentNode);
        				if(findEdge.first){
        					minValue+=liftedCosts[findEdge.second];  //add lifted edge cost if the edge exists
        				}
        				if(minValue<0){  //store only negative values
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





//    template<class LDP_STRUCT>
//    inline void ldp_single_node_cut_factor<LDP_STRUCT>::updateValues(size_t liftedEdgeID,double update){
//    	//TODO special case if number of lifted edges is low, or only one layer or whatever
//    	if(liftedEdgeID>=numberOfLiftedEdges){//update all lifted edges in the highest layer
//    		size_t lastVertex=liftedGraph.vertexFromVertex(nodeID,numberOfLiftedEdges-1);
//    		//valuesStructure[lastVertex]+=update;
//
//    		size_t highestTime=ldpStructure.getGroupIndex(lastVertex); //In both cases, maximal time for updates
//    		size_t currentLeIndex=numberOfLiftedEdges-2;
//    		size_t currentVertex=liftedGraph.vertexFromVertex(nodeID,currentLeIndex);
//    		size_t currentTime=ldpStructure.getGroupIndex(currentVertex);
//
//    		//std::unordered_set<size_t> activeVertices; //maybe only store in the
//    		//additional to valuesStructure for reconstruction of solution, maybe should be private
//
//    		while(currentTime==highestTime){
//    			if(liftedCosts[currentLeIndex]<0){
//    				valuesStructure[currentVertex]=liftedCosts[currentLeIndex]; //maybe without the plus
//    			}
//    			currentLeIndex--;
//    			if(currentLeIndex>=0){
//    				currentVertex=liftedGraph.vertexFromVertex(nodeID,currentLeIndex);
//    				currentTime=ldpStructure.getGroupIndex(currentVertex);
//    			}
//    			else{//Maybe, this was the minimal layer in the same time
//    				break;
//    			}
//    		}
//    		while(currentLeIndex>=0){ //maybe the same code for the else branch
//    			currentVertex=baseGraph.vertexFromVertex(nodeID,currentLeIndex);
//    			double minForVertex=0;
//    			size_t minIndex=ldpStructure.getTerminalNode();
//    			for (int i = 0; i < baseGraph.numberOfEdgesFromVertex(currentVertex); ++i) {
//    				size_t forwardVertex=baseGraph.vertexFromVertex(currentVertex,i);
//    				if(ldpStructure.getGroupIndex(forwardVertex)<=highestTime){//automatically reachable from nodeID
//    					//TODO if there is no lifted edge from nodeID to forward vertex, start bfs from the vertex, do not store the value
//    					//BFS: go through vertices that are not connected with nodeID via a lifted edge, maybe store as additional connections
//    					//TODO the same for cases when lifted edge cost is positive?
//    					if(valuesStructure.count(forwardVertex)>0){
//    						if(valuesStructure[forwardVertex]<minForVertex){
//    							minForVertex=valuesStructure[forwardVertex];
//    							minIndex=forwardVertex;
//    						}
//    					}
//    				}
//    			}
//    			if(minForVertex+liftedCosts[currentLeIndex]<0){ //Maybe remove this constraint for completeness? No, they will never be selected
//    				valuesStructure[currentVertex]=minForVertex+liftedCosts[currentLeIndex];
//    				indexStructure[currentVertex]=minIndex;
//    			}
//    			currentLeIndex--;
//    			if(currentLeIndex>=0){
//    				currentVertex=liftedGraph.vertexFromVertex(nodeID,currentLeIndex);
//    				currentTime=ldpStructure.getGroupIndex(currentVertex);
//    			}
//    			else{//Maybe, this was the minimal layer in the same time
//    				break;
//    			}
//    		}
//    		double minValue=0;
//    		for (int i = 0; i < numberOfEdges; ++i) {
//				size_t vertex=baseGraph.vertexFromVertex(nodeID,i);  //TODO implement several special methods calling "from" or "to" vertex based on task type
//				size_t edge=baseGraph.edgeFromVertex(nodeID,i);
//				if(valuesStructure.count(vertex)>0){
//					double value=valuesStructure[vertex]+baseCosts[edge];
//					if(value<minValue){
//						primal_=i;
//						minValue=value;
//					}
//				}
//				else{
//					double liftedECost=0;
//					auto fe=liftedGraph.findEdge(nodeID,vertex);
//					if(fe.first){
//						liftedECost=liftedCosts[fe.second];
//					}
//					double value=liftedECost+baseCosts[edge];
//					if(value<minValue){
//						primal_=i;
//						minValue=value;
//					}
//				}
//			}
//    		primalValue_=minValue;
//    	}
//    	else{//update just one edge and skip the rest in the same layer. TODO merge with the if block, most of the code is the same. Difference in the init layer
//    		size_t lastVertex=liftedGraph.vertexFromVertex(nodeID,liftedEdgeID);
//    		size_t highestTime=ldpStructure.getGroupIndex(lastVertex);
//    		valuesStructure[lastVertex]+=update;
//
//    		size_t currentLeIndex=liftedEdgeID-1;
//    		size_t currentVertex=liftedGraph.vertexFromVertex(nodeID,currentLeIndex);
//    		size_t currentTime=ldpStructure.getGroupIndex(currentVertex);
//    		while(currentTime==highestTime){
//    			currentLeIndex--;
//    			if(currentLeIndex>=0){
//    				currentVertex=liftedGraph.vertexFromVertex(nodeID,currentLeIndex);
//    				currentTime=ldpStructure.getGroupIndex(currentVertex);
//    			}
//    			else{//Maybe, this was the minimal layer in the same time
//    				break;
//    			}
//    		}
//
//    	}
//
//
//
//    }






    class ldp_mcf_single_node_cut_message
    {
        template<typename SINGLE_NODE_CUT_FACTOR>
            void RepamLeft(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
            {
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
