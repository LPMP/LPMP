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
        assert(listOfCosts.size()==listOfVertices.size());
        numberOfEdges=listOfCosts.size();
        primalSolution=std::vector<char>(numberOfEdges);
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


    void updateEdgeCost(const size_t& edgeIndex);

    double getMinMarginal(const size_t& edgeIndex) const;




private:
    double ldp_path_factor::minimize(size_t edgeIndex,bool edgeLabel)const;

    const std::vector<size_t> listOfVertices;
    std::vector<double> listOfCosts;
    const std::vector<char> isLifted;
    size_t numberOfEdges;
    std::vector<char> primalSolution;

};
double ldp_path_factor::getMinMarginal(const size_t &edgeIndex) const{
   double restrictOne=minimize(edgeIndex,1);
   double restrictZero=minimize(edgeIndex,0);
   return restrictOne-restrictZero;
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
    for (size_t i=0;i<numberOfEdges-1;i++) {
        value+=listOfCosts[i]*primalSolution[i];
    }
}

double ldp_path_factor::minimize(size_t edgeIndex,bool edgeLabel)const{
    double optimalValue=0;
    size_t indexOfWeakestNegative=numberOfEdges;
    size_t numberOfPositive=0;
    size_t indexOfPositive=numberOfEdges;
    if(edgeIndex<numberOfEdges){
        if(edgeLabel){
            optimalValue=listOfCosts[edgeIndex];
        }
        else{
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
        if(edgeIndex<numberOfEdges&&!edgeLabel){  //must cut the edge
            return optimalValue-listOfCosts[indexOfWeakestNegative];
        }
        else{ //can choose between cut one more edge or join all
            double allActive=optimalValue+listOfCosts[indexOfPositive];
            double cut=optimalValue-listOfCosts[indexOfWeakestNegative];
            return std::min(cut,allActive);
        }
    }
    else{  //no contradiction
        return optimalValue;
    }


}


double ldp_path_factor::LowerBound()const{
  double optValue=minimize(numberOfEdges,0);
  return optValue;

}

}


#endif // LDP_PATH_FACTOR_HXX
