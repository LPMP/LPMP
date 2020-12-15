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
        numberOfVerticesComplete=pInstance->getNumberOfVertices();  //without s and t

    }
    ldp_min_marginals_extractor():
    pInstance(nullptr),
    p_single_node_cut_factors_(nullptr)
    {
        numberOfVerticesComplete=0;
    }

    std::vector<std::map<size_t,double>>& getBaseEdgesMinMarginals(){
       return baseEdgesWithCosts;
    }

    std::vector<std::map<size_t,double>>& getLiftedEdgesMinMarginals(){
       return liftedEdgesWithCosts;
    }



void initMinMarginals();
void initMinMarginalsLiftedFirst();

void clearMinMarginals(){
    baseEdgesWithCosts.clear();
    liftedEdgesWithCosts.clear();
}


private:
const lifted_disjoint_paths::LdpInstance * pInstance;
std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>* p_single_node_cut_factors_;
size_t numberOfVerticesComplete;
std::vector<std::map<size_t,double>> baseEdgesWithCosts;
std::vector<std::map<size_t,double>> liftedEdgesWithCosts;
//LdpDirectedGraph baseGraph;
//LdpDirectedGraph liftedGraph;


};


template <class SINGLE_NODE_CUT_FACTOR>
inline void ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR>::initMinMarginals(){
  //  baseGraph.setAllCostToZero();
    //liftedGraph.setAllCostToZero();

    const lifted_disjoint_paths::LdpInstance &instance=*pInstance;
    assert(p_single_node_cut_factors_!=nullptr);
    std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>& single_node_cut_factors_=*p_single_node_cut_factors_;

    liftedEdgesWithCosts=std::vector<std::map<size_t,double>>(numberOfVerticesComplete-2);
    baseEdgesWithCosts=std::vector<std::map<size_t,double>>(numberOfVerticesComplete);

    //Getting the edge costs
    for (size_t i = 0; i < numberOfVerticesComplete-2; ++i) {
        // std::cout<<"node "<<i<<std::endl;

        //TODO just half of the change for base!
        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
        std::vector<double> minMarginalsIn=sncFactorIn->getAllBaseMinMarginals();
        std::vector<double> localBaseCostsIn=sncFactorIn->getBaseCosts();

        //std::cout<<"base mm size "<<minMarginalsIn.size()<<std::endl;
        for (size_t j = 0; j < minMarginalsIn.size(); ++j) {
            size_t neighborID=sncFactorIn->getBaseIDs()[j];
            //if(neighborID>=numberOfVertices-2) continue;
            minMarginalsIn[j]*=0.5;
            localBaseCostsIn.at(j)-=minMarginalsIn.at(j);
            baseEdgesWithCosts[neighborID][i]+=minMarginalsIn[j];


        }

        //  std::cout<<"base min marginals in "<<i<<std::endl;

        std::vector<double> minMarginalsLiftedIn=sncFactorIn->getAllLiftedMinMarginals(&localBaseCostsIn);
        for (size_t j = 0; j < minMarginalsLiftedIn.size(); ++j) {
            size_t neighborID=sncFactorIn->getLiftedIDs()[j];
            liftedEdgesWithCosts[neighborID][i]+=minMarginalsLiftedIn[j];

        }

          //  std::cout<<"lifted min marginals in "<<i<<std::endl;


        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();
        std::vector<double> localBaseCostsOut=sncFactorOut->getBaseCosts();
        std::vector<double> minMarginalsOut=sncFactorOut->getAllBaseMinMarginals();

        for (size_t j = 0; j < minMarginalsOut.size(); ++j) {
            size_t neighborID=sncFactorOut->getBaseIDs()[j];
            //if(neighborID>=numberOfVertices-2) continue;
            minMarginalsOut[j]*=0.5;
              localBaseCostsOut.at(j)-=minMarginalsOut.at(j);
            baseEdgesWithCosts[i][neighborID]+=minMarginalsOut[j];

        }

         std::vector<double> minMarginalsLiftedOut=sncFactorOut->getAllLiftedMinMarginals(&localBaseCostsOut);

         for (size_t j = 0; j < minMarginalsLiftedOut.size(); ++j) {
             size_t neighborID=sncFactorOut->getLiftedIDs()[j];
             liftedEdgesWithCosts[i][neighborID]+=minMarginalsLiftedOut[j];
         }

    }


}


template <class SINGLE_NODE_CUT_FACTOR>
inline void ldp_min_marginals_extractor<SINGLE_NODE_CUT_FACTOR>::initMinMarginalsLiftedFirst(){
  //  baseGraph.setAllCostToZero();
    //liftedGraph.setAllCostToZero();

    const lifted_disjoint_paths::LdpInstance &instance=*pInstance;
    assert(p_single_node_cut_factors_!=nullptr);
    std::vector<std::array<SINGLE_NODE_CUT_FACTOR*,2>>& single_node_cut_factors_=*p_single_node_cut_factors_;

    liftedEdgesWithCosts=std::vector<std::map<size_t,double>>(numberOfVerticesComplete-2);
    baseEdgesWithCosts=std::vector<std::map<size_t,double>>(numberOfVerticesComplete);

    //Getting the edge costs
    for (size_t i = 0; i < numberOfVerticesComplete-2; ++i) {
        // std::cout<<"node "<<i<<std::endl;

        //TODO just half of the change for base!
        const auto* sncFactorIn=single_node_cut_factors_[i][0]->get_factor();
        std::vector<double> liftedMinMarginalsIn=sncFactorIn->getAllLiftedMinMarginals();
        std::vector<double> localLiftedCostsIn=sncFactorIn->getLiftedCosts();

        assert(localLiftedCostsIn.size()==liftedMinMarginalsIn.size());


        //std::vector<double> minMarginalsLiftedIn=sncFactorIn->getAllLiftedMinMarginals(&localBaseCostsIn);
        for (size_t j = 0; j < liftedMinMarginalsIn.size(); ++j) {
            liftedMinMarginalsIn[j]*=0.5;
            localLiftedCostsIn[j]-=liftedMinMarginalsIn[j];
            size_t neighborID=sncFactorIn->getLiftedIDs()[j];
            liftedEdgesWithCosts[neighborID][i]+=liftedMinMarginalsIn[j];

        }

        const std::vector<double>& baseCostsIn=sncFactorIn->getBaseCosts();
        std::vector<double> baseMinMarginalsIn=sncFactorIn->getAllBaseMinMarginals(&baseCostsIn,&localLiftedCostsIn);

        //std::cout<<"base mm size "<<minMarginalsIn.size()<<std::endl;
        for (size_t j = 0; j < baseMinMarginalsIn.size(); ++j) {
            size_t neighborID=sncFactorIn->getBaseIDs()[j];
            baseEdgesWithCosts[neighborID][i]+=baseMinMarginalsIn[j];
        }



        const auto* sncFactorOut=single_node_cut_factors_[i][1]->get_factor();
        std::vector<double> liftedMinMarginalsOut=sncFactorOut->getAllLiftedMinMarginals();
        std::vector<double> localLiftedCostsOut=sncFactorOut->getLiftedCosts();

        for (size_t j = 0; j < liftedMinMarginalsOut.size(); ++j) {
            liftedMinMarginalsOut[j]*=0.5;
            localLiftedCostsOut[j]-=liftedMinMarginalsOut[j];
            size_t neighborID=sncFactorOut->getLiftedIDs()[j];
            liftedEdgesWithCosts[i][neighborID]+=liftedMinMarginalsOut[j];
        }

        const std::vector<double>& baseCostsOut=sncFactorOut->getBaseCosts();
        std::vector<double> baseMinMarginalsOut=sncFactorOut->getAllBaseMinMarginals(&baseCostsOut,&localLiftedCostsOut);



        for (size_t j = 0; j < baseMinMarginalsOut.size(); ++j) {
            size_t neighborID=sncFactorOut->getBaseIDs()[j];
            baseEdgesWithCosts[i][neighborID]+=baseMinMarginalsOut[j];

        }


    }


}




}


#endif // LDP_MIN_MARGINALS_EXTRACTOR_HXX
