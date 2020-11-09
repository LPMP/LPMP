#ifndef LDP_PATH_FACTOR_HXX
#define LDP_PATH_FACTOR_HXX

#include<cstdlib>
#include<vector>
#include <config.hxx>

namespace LPMP {

class ldp_path_factor{

public:
    ldp_path_factor(const std::vector<size_t>& _listOfVertices,const std::vector<double>& _listOfCosts,const std::vector<char>& _isLifted):
      listOfVertices(_listOfVertices),
      listOfCosts(_listOfCosts),
      isLifted(_isLifted)
    {
        assert(listOfCosts.size()==listOfVertices.size()&&listOfCosts.size()==isLifted.size());
        numberOfEdges=listOfCosts.size();
        primalSolution=std::vector<char>(numberOfEdges);
        assert(listOfCosts.size()>=2);
    }



    double LowerBound() const;

    double EvaluatePrimal() const;

    void setPrimal(const std::vector<size_t>& primalDescendants, const std::vector<size_t> &vertexLabels);

    const std::vector<char>& getPrimal(){
        return primalSolution;
    }

    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

    auto export_variables() { return std::tie(); }

    void init_primal();

   // std::array<double,3> getAllMinMarginals();


    void updateEdgeCost(const size_t& edgeIndex,const double& value);

    double getMinMarginal(const size_t& edgeIndex) const;


    const std::vector<size_t>& getListOfVertices()const{
        return listOfVertices;
    }

    const std::vector<double>& getCosts()const{
        return listOfCosts;
    }

    const std::vector<char>& getLiftedInfo()const{
        return isLifted;
    }

    const size_t& getNumberOfEdges()const{
        return numberOfEdges;
    }


private:
    double minimize(size_t edgeIndex,bool edgeLabel)const;

    const std::vector<size_t> listOfVertices;
    std::vector<double> listOfCosts;
    const std::vector<char> isLifted;
    size_t numberOfEdges;
    std::vector<char> primalSolution;

};


double ldp_path_factor::getMinMarginal(const size_t &edgeIndex) const{
      assert(listOfCosts.size()==listOfVertices.size()&&listOfCosts.size()==isLifted.size());
      assert(listOfCosts.size()==numberOfEdges);
      assert(numberOfEdges>=2);
    assert(edgeIndex<numberOfEdges);
    std::cout<<" printing path factor "<<std::endl;
    for(size_t i=0;i<numberOfEdges;i++){
        std::cout<<"vertex "<<listOfVertices[i]<<std::endl;
        //std::cout<<"vertex "<<listOfVertices[i]<<std::endl;
    }

   double restrictOne=minimize(edgeIndex,1);
   std::cout<<"restrict one "<<restrictOne<<std::endl;
   double restrictZero=minimize(edgeIndex,0);
   std::cout<<"restrict zero "<<restrictZero<<std::endl;
   return restrictOne-restrictZero;
}

void ldp_path_factor::init_primal(){
    std::fill(primalSolution.begin(),primalSolution.end(),0);
}

void ldp_path_factor::updateEdgeCost(const size_t& edgeIndex,const double& value){
    assert(edgeIndex<listOfCosts.size());
    std::cout<<"update edge cost, index "<<edgeIndex<<", value "<<value<<std::endl;
    listOfCosts[edgeIndex]+=value;
}

void ldp_path_factor::setPrimal(const std::vector<size_t>& primalDescendants, const std::vector<size_t> &vertexLabels){

    size_t numberOfZeros=0;
    for (size_t i=0;i<numberOfEdges-1;i++) {

        size_t vertex1=listOfVertices[i];
        size_t vertex2=listOfVertices[i+1];

        if(isLifted[i]){
            if(vertexLabels[vertex1]==vertexLabels[vertex2]){
                primalSolution[i]=1;
            }
            else{
                primalSolution[i]=0;
                numberOfZeros++;
            }
        }
        else{
            if(primalDescendants[vertex1]==vertex2){
                primalSolution[i]=1;
            }
            else {
                primalSolution[i]=0;
                numberOfZeros++;
            }

        }
    }
    size_t vertex2=listOfVertices.back();
    size_t vertex1=listOfVertices[0];
    primalSolution.back()=(vertexLabels[vertex1]==vertexLabels[vertex2]);
    if(primalSolution.back()==0) numberOfZeros++;
    assert(numberOfZeros!=1);

}

double ldp_path_factor::EvaluatePrimal() const{
    double value=0;
    for (size_t i=0;i<numberOfEdges;i++) {
        value+=listOfCosts[i]*primalSolution[i];
    }
    return value;
}

double ldp_path_factor::minimize(size_t edgeIndex,bool edgeLabel)const{
    double optimalValue=0;
    size_t indexOfWeakestNegative=numberOfEdges;
    size_t numberOfPositive=0;
    size_t indexOfPositive=numberOfEdges;
    std::cout<<"calling minimize, printing edge costs "<<std::endl;
    for(size_t i=0;i<numberOfEdges;i++){
        double c=listOfCosts[i];
        size_t vertex=listOfVertices[i];
        std::cout<<listOfVertices[i]<<": "<<c<<std::endl;

    }
    if(edgeIndex<numberOfEdges){
//        std::cout<<"printing edge costs "<<std::endl;
//        for(auto& c:listOfCosts){
//            std::cout<<c<<std::endl;
//        }

        if(edgeLabel){
            std::cout<<"minimize with restrict one "<<edgeIndex<<std::endl;
            optimalValue=listOfCosts[edgeIndex];
        }
        else{
            std::cout<<"minimize with restrict zero"<<edgeIndex<<std::endl;
            numberOfPositive=1;
            indexOfPositive=edgeIndex;
        }
    }


    for(size_t i=0;i<numberOfEdges;i++){
        if(i==edgeIndex) continue;
        if(listOfCosts[i]<0){
            optimalValue+=listOfCosts[i];
            if(indexOfWeakestNegative==numberOfEdges){
                indexOfWeakestNegative=i;
            }
            else if(listOfCosts[i]>listOfCosts[indexOfWeakestNegative]){
                indexOfWeakestNegative=i;
            }
        }
        else{
            numberOfPositive++;
            if(numberOfPositive==1){
                indexOfPositive=i;
            }
        }

    }
    if(numberOfPositive==1){
        assert(indexOfWeakestNegative<numberOfEdges);
        assert(indexOfPositive<numberOfEdges);
        if(edgeIndex<numberOfEdges&&!edgeLabel){  //must cut the edge

            double value=optimalValue-listOfCosts[indexOfWeakestNegative];
            std::cout<<"restrict zero contradiction, return  "<<value<<std::endl;
            return value;
        }
        else{ //can choose between cut one more edge or join all

            double allActive=optimalValue+listOfCosts[indexOfPositive];
            if(edgeIndex<numberOfEdges)std::cout<<"all active cost "<<allActive<<std::endl;
            double cut=optimalValue-listOfCosts[indexOfWeakestNegative];
            if(edgeIndex<numberOfEdges)std::cout<<"do cut cost "<<cut<<std::endl;
            return std::min(cut,allActive);
        }
    }
    else{  //no contradiction
        if(edgeIndex<numberOfEdges) std::cout<<"no cotradiction "<<std::endl;
        return optimalValue;
    }


}


double ldp_path_factor::LowerBound()const{
  double optValue=minimize(numberOfEdges,0);
  return optValue;

}


class ldp_snc_path_message
{
public:
    ldp_snc_path_message(const std::vector<size_t>& _edgeIndexInPath,  //Mostly one, two for the first and the last path vertices
         const std::vector<size_t>& _vertexIndexInSnc,
         const std::vector<char>& _isLifted)
        {
        edgeIndexInPath=_edgeIndexInPath;  //Mostly one, two for the first and the last path vertices
         vertexIndexInSnc=_vertexIndexInSnc;
        isLifted=_isLifted;
            dimension=_edgeIndexInPath.size();
            std::cout<<"dimension in message "<<dimension<<std::endl;

            std::cout<<"dimension of edge indices "<<edgeIndexInPath.size()<<std::endl;
           assert(dimension==1||dimension==2);
          assert(vertexIndexInSnc.size()==dimension);
            assert(isLifted.size()==dimension);
            assert(edgeIndexInPath.size()==dimension);

        }

    template<typename PATH_FACTOR>
    void RepamLeft(PATH_FACTOR& l, const double msg, const std::size_t msg_dim) const
    {
      std::cout<<"repam left path "<<std::endl;
      double before=l.LowerBound();
      std::cout<<"before lb "<<before<<std::endl;
      assert(dimension==1||dimension==2);
      assert(msg_dim<dimension);
     assert(vertexIndexInSnc.size()==dimension);
       assert(isLifted.size()==dimension);
       assert(edgeIndexInPath.size()==dimension);
       l.updateEdgeCost(edgeIndexInPath[msg_dim],msg);
       double after=l.LowerBound();
       std::cout<<"after lb "<<after<<std::endl;
     }

    template<typename SINGLE_NODE_CUT_FACTOR>
    void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
    {
        std::cout<<"snc repam right "<<std::endl;
        double before=r.LowerBound();
        std::cout<<"before lb "<<before<<std::endl;
        assert(dimension==1||dimension==2);
        assert(msg_dim<dimension);
       assert(vertexIndexInSnc.size()==dimension);
         assert(isLifted.size()==dimension);
         assert(edgeIndexInPath.size()==dimension);
              std::cout<<"repam left path "<<std::endl;
        r.updateEdgeCost(msg,vertexIndexInSnc[msg_dim],isLifted[msg_dim]);
        double after=r.LowerBound();
        std::cout<<"abter lb "<<after<<std::endl;
    }

    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
    void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
    {
        assert(dimension==1||dimension==2);
       assert(vertexIndexInSnc.size()==dimension);
         assert(isLifted.size()==dimension);
         assert(edgeIndexInPath.size()==dimension);
           std::cout<<"send messages to left path "<<std::endl;
              assert(isLifted.size()==dimension);
        assert(dimension==1||dimension==2);
        double delta=0;
        for (size_t i=0;i<dimension;i++) {
            if(isLifted.at(i)){
                         std::cout<<"lifted "<<std::endl;

                delta = r.getOneLiftedMinMarginal(vertexIndexInSnc[i]);
                // std::cout<<"delta obtained "<<std::endl;
            }
            else{
                std::cout<<"base "<<std::endl;
                delta = r.getOneBaseEdgeMinMarginal(vertexIndexInSnc[i]);

            }


            msg[i] -= omega * delta;
        }

    }

    template<typename PATH_FACTOR, typename MSG>
    void send_message_to_right(const PATH_FACTOR& l, MSG& msg, const double omega)
    {        std::cout<<"send messages to right path "<<std::endl;
        assert(dimension==1||dimension==2);
             assert(dimension==1||dimension==2);
            assert(vertexIndexInSnc.size()==dimension);
              assert(isLifted.size()==dimension);
              assert(edgeIndexInPath.size()==dimension);


        double delta;

        for (size_t i=0;i<dimension;i++) {
            delta=l.getMinMarginal(edgeIndexInPath[i]);
            std::cout<<"sending "<<delta<<std::endl;
            msg[i] -= omega * delta;
        }

    }


    template<typename SINGLE_NODE_CUT_FACTOR,typename CUT_FACTOR>
    bool check_primal_consistency(const CUT_FACTOR& l, const SINGLE_NODE_CUT_FACTOR& r) const
    {
        return true;

    }

private:
    std::vector<size_t> edgeIndexInPath;  //Mostly one, two for the first and the last path vertices
    size_t dimension;
    std::vector<size_t> vertexIndexInSnc;
    std::vector<char> isLifted;

};




}


#endif // LDP_PATH_FACTOR_HXX
