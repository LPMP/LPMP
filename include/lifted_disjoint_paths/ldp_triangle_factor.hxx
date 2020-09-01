/*
 * ldp_triangle_factor.hxx
 *
 *  Created on: Apr 6, 2020
 *      Author: fuksova
 */

#ifndef INCLUDE_LIFTED_DISJOINT_PATHS_LDP_TRIANGLE_FACTOR_HXX_
#define INCLUDE_LIFTED_DISJOINT_PATHS_LDP_TRIANGLE_FACTOR_HXX_

#include <cstdlib>
#include <string>         // std::string
#include <bitset>
#include <vector>
#include<array>
#include <config.hxx>
#include <assert.h>
#include <unordered_map>
namespace LPMP {

class ldp_triangle_factor
{
public:

	//By default, all edges are lifted. However vu or uw can be base too.
    ldp_triangle_factor(size_t v,size_t u,size_t w,const std::array<double,3>& costs,bool vuBase=0,bool uwBase=0): //TODO remember vertex indices instead of edge indices?
		vInd(v),
		uInd(u),
        wInd(w),
        v1v2Base(vuBase),
        v2v3Base(uwBase),
        edgeCosts(costs)
{
//        assert(costs.size()==3);

//        for (int i = 0; i < 3; ++i) {
//            edgeCosts[i]=costs[i];
//        }
		//Feasible labelings of edge (vu,uw,vw)
		labelings={std::bitset<3>("000"),std::bitset<3>("111"),std::bitset<3>("100"),std::bitset<3>("010"),std::bitset<3>("001")};
		if(vuBase){
			labelings.push_back(std::bitset<3>("011"));
		}
		if(uwBase){
			labelings.push_back(std::bitset<3>("101"));
		}
		primal_=0;
}

    ldp_triangle_factor(const ldp_triangle_factor& trFactor): //TODO remember vertex indices instead of edge indices?
        vInd(trFactor.vInd),
        uInd(trFactor.uInd),
        wInd(trFactor.wInd),
        v1v2Base(trFactor.v1v2Base),
        v2v3Base(trFactor.v2v3Base),
        edgeCosts(trFactor.edgeCosts),
        labelings(trFactor.labelings),
        primal_(trFactor.primal_)
{

}


    const std::array<double,3> getEdgeCosts()const {
        return edgeCosts;
    }

    const size_t & getV1()const {
        return vInd;
    }


    const size_t & getV2()const {
        return uInd;
    }


    const size_t & getV3()const {
        return wInd;
    }

    bool isV1V2Base()const{
        return v1v2Base;
    }

    bool isV2V3Base()const{
        return v2v3Base;
    }

    void print() const{
        std::cout<<vInd<<","<<uInd<<","<<wInd<<" ";
    }
	double LowerBound() const{
		double minValue=0;
		short minIndex=0;
		for (short i = 0; i < labelings.size(); ++i) {
            const std::bitset<3>& label=labelings.at(i);
			double value=label[0]*edgeCosts[0]+label[1]*edgeCosts[1]+label[2]*edgeCosts[2];
			if(value<minValue){
				minValue=value;
				minIndex=i;
			}
		}
		return minValue;
	}

	double EvaluatePrimal() const{
        const std::bitset<3>& label=labelings.at(primal_);
		double value=label[0]*edgeCosts[0]+label[1]*edgeCosts[1]+label[2]*edgeCosts[2];
		return value;
	}

    void setPrimal(const std::bitset<3>& primalSolution){
        bool primalSet=false;
        for (short i = 0; i < labelings.size(); ++i) {
            const std::bitset<3>& label=labelings.at(i);
            if(label==primalSolution){
                primal_=i;
                primalSet=true;
                break;
            }
        }
        assert(primalSet);
    }

    const std::bitset<3>& getPrimal()const{
        return labelings.at(primal_);
    }

    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }
    // auto export_variables() { return std::tie(); }

    auto export_variables() { return std::tie(); }

	void init_primal() { primal_ = 0; }

    double delta(short edgeId,const std::array<double,3>& localEdgeCosts)const{  //difference: label[edgeId]=1-label[edgeId]=0
		assert(edgeId<=2);
		double min1=std::numeric_limits<double>::infinity();
		double min0=std::numeric_limits<double>::infinity();
		for (int i = 0; i < labelings.size(); ++i) {
            double value=labelings[i][0]*localEdgeCosts[0]+labelings[i][1]*localEdgeCosts[1]+labelings[i][2]*localEdgeCosts[2];
			if(labelings[i][edgeId]){
				min1=std::min(min1,value);
			}
			else{
				min0=std::min(min0,value);
			}
		}
		return min1-min0;
	}

    std::array<double,3> getAllMinMarginals() const{
        std::array<double,3> localEdgeCosts=edgeCosts;
        std::array<double,3> minMarginals;
        double d=delta(0,localEdgeCosts);
        localEdgeCosts[0]-=d;
        minMarginals[0]=d;

        d=delta(1,localEdgeCosts);
        localEdgeCosts[1]-=d;
        minMarginals[1]=d;

        d=delta(2,localEdgeCosts);
        localEdgeCosts[2]-=d;
        minMarginals[2]=d;

        return minMarginals;
	}

    double getOneMinMarginal(size_t edgeId) const{
       assert(edgeId<=2);
        double d=delta(edgeId,edgeCosts);
        return d;
    }

    void updateCost(size_t edgeId,double update) {  //update cost of one edge, assumed local edge ID: indices 0-2
        if(edgeId>2){
            std::cout<<"error in triangle update cost, edge id "<<edgeId<<std::endl;
        }
		assert(edgeId<=2);
		edgeCosts[edgeId]+=update;
	}

//    void updateCost(size_t v1, size_t v2,double update){ //assumes directed edge vertices
//        size_t e=getLocalEdgeID(v1,v2);
//		updateCost(e,update);
//	}

	size_t getLocalEdgeID(size_t v1, size_t v2){ //assumes directed edge vertices. If edge not present in the factor, returns 3.
		if(v1==vInd){
			if(v2==uInd){
				return 0;
			}
			else if(v2==wInd){
				return 2;
			}
		}
		else if(v1==uInd){
			if(v2==wInd){
				return 1;
			}
		}
        assert(false);
	}

private:
	std::size_t primal_; // index of feasible labeling
	size_t vInd,uInd,wInd;
    //double edgeCosts[3]; //vu, uw, vw
    std::array<double,3> edgeCosts;
	std::vector<std::bitset<3>> labelings;
    const bool v1v2Base;
    const bool v2v3Base;
};





class ldp_snc_triangle_message
{
public:
    ldp_snc_triangle_message(const std::vector<size_t>& _indicesInTriangle, const std::vector<size_t>& _verticesInSnc,const std::vector<bool> & _isLifted)
    : verticesInSnc(_verticesInSnc),
      indicesInTriangle(_indicesInTriangle),
      isLifted(_isLifted)
    {
        assert(verticesInSnc.size()==isLifted.size()&&indicesInTriangle.size()==isLifted.size());
        assert(verticesInSnc.size()<=2);
        for (int i = 0; i < verticesInSnc.size(); ++i) {
            assert(indicesInTriangle[i]<=2);
        }
    }

    template<typename TRIANGLE_FACTOR>
    void RepamLeft(TRIANGLE_FACTOR& l, const double msg, const std::size_t msg_dim) const
    {
       if(debug()) std::cout<<"triangle repam left ";
//       l.print();
//       printIndices();

      // std::cout<<std::endl;
        assert(msg_dim <=indicesInTriangle.size());
       //std::cout<<"index to update "<<indicesInTriangle.at(msg_dim);
        l.updateCost(indicesInTriangle.at(msg_dim),msg);
     }

    template<typename SINGLE_NODE_CUT_FACTOR>
    void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
    {
        if(debug()) std::cout<<"triangle repam right "<<std::endl;
        assert(msg_dim == 0||msg_dim==1);
        r.updateCostSimple(msg,verticesInSnc.at(msg_dim),isLifted.at(msg_dim));
    }

    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
    void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
    {
        if(debug()) std::cout<<"triangle send one to left ";
        //printIndices();
        //for (size_t i=0;i<1;i++) {
        for (size_t i=0;i<isLifted.size();i++) {
            double delta=0;
            if(isLifted[i]){
               // std::cout<<"lifted "<<std::endl;
                delta = r.oneLiftedMinMarginal(verticesInSnc.at(i));
            }
            else{
              //  std::cout<<"base "<<std::endl;
                delta = r.getOneBaseEdgeMinMarginal(verticesInSnc.at(i));
            }
            msg[i] -= omega * delta;
        }
      //  std::cout<<" done "<<std::endl;
    }

    template<typename TRIANGLE_FACTOR, typename MSG>
    void send_message_to_right(const TRIANGLE_FACTOR& l, MSG& msg, const double omega)
    {
        if(debug()) std::cout<<"triangle send one to right"<<std::endl;
        //for (size_t i=0;i<1;i++) {
        for (size_t i=0;i<indicesInTriangle.size();i++) {
            const double delta = l.getOneMinMarginal(indicesInTriangle.at(i));
            msg[i] -= omega * delta;
        }
    }



    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG_ARRAY>
    static void SendMessagesToLeft(const SINGLE_NODE_CUT_FACTOR& r, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
    {
        if(debug()) std::cout<<"triangle send all to left"<<std::endl;

        //TODO take only one half from the base min marginals!
        std::vector<double> msg_vec_base = r.getAllBaseMinMarginals();
        std::vector<double> localBaseCost=r.getBaseCosts();
        for (size_t i = 0; i < localBaseCost.size(); ++i) {
            msg_vec_base[i]=0.5*msg_vec_base[i];
            localBaseCost[i]-=msg_vec_base[i];
        }
        const std::vector<double> msg_vec_lifted = r.getAllLiftedMinMarginals(&localBaseCost);

        for(auto it=msg_begin; it!=msg_end; ++it)
        {
            auto& msg = (*it).GetMessageOp();
            for (int i = 0; i <msg.verticesInSnc.size(); ++i) {
                const size_t vertex = msg.verticesInSnc.at(i);
                bool lifted=msg.isLifted.at(i);
                if(lifted){
                    (*it)[i] -=  omega * msg_vec_lifted.at(vertex);
                }
                else{
                    (*it)[i] -= omega * msg_vec_base.at(vertex);
                }
//                if(lifted){
//                    (*it)[i] -= r.getTriangleCoeffForLiftedEdge(vertex)* omega * msg_vec_lifted.at(vertex);
//                }
//                else{
//                    (*it)[i] -= r.getTriangleCoeffForBaseEdge(vertex)*omega * msg_vec_base.at(vertex);
//                }
            }
        }
    }

    template<typename TRIANGLE_FACTOR, typename MSG_ARRAY>
    static void SendMessagesToRight(const TRIANGLE_FACTOR& l, MSG_ARRAY msg_begin, MSG_ARRAY msg_end, const double omega)
    {

        if(debug()) std::cout<<"triangle send all to right"<<std::endl;
        const std::array<double,3> msg_vec = l.getAllMinMarginals();

        for(auto it=msg_begin; it!=msg_end; ++it)
        {
            auto& msg = (*it).GetMessageOp();
            for (int i = 0; i < msg.indicesInTriangle.size(); ++i) {
                const size_t indexInTr = msg.indicesInTriangle.at(i);
                //(*it)[i]-= 0.5*omega * msg_vec.at(indexInTr);
                (*it)[i]-= omega * msg_vec.at(indexInTr);
            }

        }

    }



    template<typename SINGLE_NODE_CUT_FACTOR,typename TRIANGLE_FACTOR>
    bool check_primal_consistency(const TRIANGLE_FACTOR& l, const SINGLE_NODE_CUT_FACTOR& r) const
    {
        return true;
//        if(debug()) std::cout<<"triangle check primal"<<std::endl;
//        bool isConsistent=true;
//        for (int i = 0; i < verticesInSnc.size()&&isConsistent; ++i) {
//        //for (int i = 0; i < 1&&isConsistent; ++i) {
//            const bool triangleActive = l.getPrimal().at(indicesInTriangle.at(i));
//            bool edgeActive=false;
//            if(isLifted[i]){
//                edgeActive=r.isActiveInPrimalLifted(verticesInSnc.at(i));
//            }
//            else{
//                edgeActive=r.getPrimalBaseIndex()==verticesInSnc.at(i);
//            }
//            isConsistent=(edgeActive==triangleActive);
//        }
//        return isConsistent;
    }

private:
//    void printIndices()const{
//        std::cout<<"indices in triangle ";
//        for (int i = 0; i < indicesInTriangle.size(); ++i) {
//            std::cout<<indicesInTriangle.at(i)<<",";
//        }
//        std::cout<<std::endl;
//    }
    const std::vector<size_t> indicesInTriangle;
    const std::vector<size_t> verticesInSnc;
    const std::vector<bool> isLifted;
};




}





#endif /* INCLUDE_LIFTED_DISJOINT_PATHS_LDP_TRIANGLE_FACTOR_HXX_ */
