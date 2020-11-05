#ifndef LDP_MIN_MARGINALS_EXTRACTOR_HXX
#define LDP_MIN_MARGINALS_EXTRACTOR_HXX

#include"lifted_disjoint_paths/ldp_instance.hxx"


namespace LPMP {

template <class SINGLE_NODE_CUT_FACTOR>
class ldp_min_marginals_extractor {

public:
    ldp_min_marginals_extractor(std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>* _p_sncFactorContainer, const lifted_disjoint_paths::LdpInstance * _pInstance)//:
  // pInstance(_pInstance),
   // p_single_node_cut_factors_(_p_sncFactorContainer)
    {
        p_single_node_cut_factors_=_p_sncFactorContainer;
        pInstance=_pInstance;
       // baseGraph=pInstance->getMyGraph();
       // liftedGraph=pInstance->getMyGraphLifted();
        //directedGraph.setAllCostToZero();
        numberOfVertices=pInstance->getNumberOfVertices()-2;  //without s and t

    }
    ldp_min_marginals_extractor():
    pInstance(nullptr),
    p_single_node_cut_factors_(nullptr)
    {
        numberOfVertices=0;
    }

    std::vector<std::map<size_t,double>>& getBaseEdgesMinMarginals(){
       return baseEdgesWithCosts;
    }

    std::vector<std::map<size_t,double>>& getLiftedEdgesMinMarginals(){
       return liftedEdgesWithCosts;
    }



void initMinMarginals();


private:
const lifted_disjoint_paths::LdpInstance * pInstance;
std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>* p_single_node_cut_factors_;
size_t numberOfVertices;
std::vector<std::map<size_t,double>> baseEdgesWithCosts;
std::vector<std::map<size_t,double>> liftedEdgesWithCosts;
//LdpDirectedGraph baseGraph;
//LdpDirectedGraph liftedGraph;


};


template <class SINGLE_NODE_CUT_FACTOR>
void ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR>::initMinMarginals(){
  //  baseGraph.setAllCostToZero();
    //liftedGraph.setAllCostToZero();

    const lifted_disjoint_paths::LdpInstance &instance=*pInstance;
    std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>& single_node_cut_factors_=*p_single_node_cut_factors_;

    std::map<size_t,std::set<size_t>> baseEdgeUsed;
    std::map<size_t,std::set<size_t>> liftedEdgeUsed;


   // std::map<size_t,std::map<size_t,double>> baseEdgesWithCosts;  //TODO: can be vector of maps
    std::vector<std::vector<size_t>> baseEdgesIn(numberOfVertices);  //can be two dim array
    std::vector<std::vector<size_t>> baseEdgesOut(numberOfVertices);

    std::vector<std::vector<size_t>> liftedEdgesIn(numberOfVertices);
    std::vector<std::vector<size_t>> liftedEdgesOut(numberOfVertices);
    //std::multimap<double,ldp_triangle_factor> candidateFactors;

    std::vector<std::tuple<double,size_t,size_t>> edgesToSort;


    andres::graph::Digraph<> connectivityGraph(numberOfVertices);
    std::vector<std::set<size_t>> descendants(numberOfVertices);
    std::vector<std::set<size_t>> predecessors(numberOfVertices);

    liftedEdgesWithCosts=std::vector<std::map<size_t,double>>(numberOfVertices);
    baseEdgesWithCosts=std::vector<std::map<size_t,double>>(numberOfVertices);

    //Getting the edge costs
    for (size_t i = 0; i < numberOfVertices; ++i) {
        // std::cout<<"node "<<i<<std::endl;

        //TODO just half of the change for base!
        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
        std::vector<double> minMarginalsIn=sncFactorIn->getAllBaseMinMarginals();
        std::vector<double> localBaseCostsIn=sncFactorIn->getBaseCosts();

        //std::cout<<"base mm size "<<minMarginalsIn.size()<<std::endl;
        for (size_t j = 0; j < minMarginalsIn.size(); ++j) {
            size_t neighborID=sncFactorIn->getBaseIDs()[j];
            if(neighborID>=numberOfVertices) continue;
           // std::cout<<"j "<<j<<std::endl;
            minMarginalsIn[j]*=0.5;
           // std::cout<<"1"<<std::endl;
            baseEdgesIn[i].push_back(minMarginalsIn.at(j));
           // std::cout<<"1"<<std::endl;
            localBaseCostsIn.at(j)-=minMarginalsIn.at(j);
            //std::cout<<"1"<<std::endl;
            baseEdgesWithCosts[neighborID][i]+=minMarginalsIn[j];
           // std::cout<<"get the cost"<<std::endl;
        }

        //  std::cout<<"base min marginals in "<<i<<std::endl;

        std::vector<double> minMarginalsLiftedIn=sncFactorIn->getAllLiftedMinMarginals(&localBaseCostsIn);
        for (size_t j = 0; j < minMarginalsLiftedIn.size(); ++j) {
           liftedEdgesIn[i].push_back(minMarginalsLiftedIn[j]);
           liftedEdgesWithCosts[sncFactorIn->getLiftedIDs()[j]][i]+=minMarginalsLiftedIn[j];
        }

          //  std::cout<<"lifted min marginals in "<<i<<std::endl;


        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();
        std::vector<double> localBaseCostsOut=sncFactorOut->getBaseCosts();
        std::vector<double> minMarginalsOut=sncFactorOut->getAllBaseMinMarginals();

        for (size_t j = 0; j < minMarginalsOut.size(); ++j) {
            size_t neighborID=sncFactorOut->getBaseIDs()[j];
            if(neighborID>=numberOfVertices) continue;
            minMarginalsOut[j]*=0.5;
            baseEdgesOut[i].push_back(minMarginalsOut.at(j));
            localBaseCostsOut.at(j)-=minMarginalsOut.at(j);
            baseEdgesWithCosts[i][neighborID]+=minMarginalsOut[j];
        }

         std::vector<double> minMarginalsLiftedOut=sncFactorOut->getAllLiftedMinMarginals(&localBaseCostsOut);

        for (size_t j = 0; j < minMarginalsLiftedOut.size(); ++j) {
            liftedEdgesOut[i].push_back(minMarginalsLiftedOut[j]);
            liftedEdgesWithCosts[i][sncFactorOut->getLiftedIDs()[j]]+=minMarginalsLiftedOut[j];
        }

     }


}



}


#endif // LDP_MIN_MARGINALS_EXTRACTOR_HXX
