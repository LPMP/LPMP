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

            //virtual size_t getNeighborBase()

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
        private:
            void updateValues(size_t liftedEdgeID=numberOfLiftedEdges);  //Highest number: update all.For base edge update, simpler procedure
            //TODO implement update for base edge cost

            std::size_t nodeID;
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
        size_t vertexToReach=ldpStructure.getTerminalNode();
        if(oneEdgeUpdate){
        	vertexToReach=liftedGraph.edgeFromVertex(nodeID,liftedEdgeID);
        }

    	std::stack<size_t> nodeStack;
           	nodeStack.push(nodeID);

        	while(!nodeStack.empty()){
        		size_t currentNode=nodeStack.top();
        		bool descClosed=true;
        		double minValue=0;
        		size_t minValueIndex=ldpStructure.getTerminalNode;

        		for (int i = 0; i < baseGraph.numberOfEdgesFromVertex(currentNode); ++i) {
        			size_t desc=baseGraph.vertexFromVertex(currentNode,i);
        			if(ldpStructure.getGroupIndex(desc)<=maxLayer){
        				if(indexStructure.count(desc)>0||(oneEdgeUpdate&&!ldpStructure.isReachable(currentNode,vertexToReach))){  //descendant closed
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
//        						size_t forwardVertex=indexStructure[vertex];
//        						if(forwardVertex!=ldpStructure.getTerminalNode){
//        							valueToAdd+=valuesStructure[forwardVertex];
//        						}
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
        				auto findEdgeBase=baseGraph.findEdge(nodeID,currentNode);
        				if(minValue<0||(findEdgeBase.first&&minValue<liftedCosts[findEdge.second])){  //store only negative values or values needed to correct solutionCosts
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
