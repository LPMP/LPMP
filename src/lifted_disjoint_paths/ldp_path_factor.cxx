#include"lifted_disjoint_paths/ldp_path_factor.hxx"

namespace LPMP {



ldp_path_factor::ldp_path_factor(const std::vector<size_t>& _listOfVertices, const std::vector<double>& _listOfCosts, const std::vector<char>& _isLifted, const lifted_disjoint_paths::LdpInstance* pInstance, bool _mustCut):
    listOfVertices(_listOfVertices),
    listOfCosts(_listOfCosts),
    isLifted(_isLifted),
    mustCut(_mustCut)
{
    assert(!mustCut);
    if(mustCut){
        assert(listOfCosts.size()+1==listOfVertices.size());
        assert(listOfCosts.size()>=2);
    }
    else{
        assert(listOfCosts.size()==listOfVertices.size());
        assert(isLifted.back()==1);
        assert(listOfCosts.size()>2);
    }

   // assert(listOfCosts.size()==listOfVertices.size()&&listOfCosts.size()==isLifted.size());
    assert(listOfCosts.size()==isLifted.size());
    numberOfEdges=listOfCosts.size();
    isStrongBase=std::vector<char>(numberOfEdges);
    primalSolution=std::vector<char>(numberOfEdges);
    optSolution=std::vector<char>(numberOfEdges);
   // assert(listOfCosts.size()>2);
    for (size_t i = 0; i < numberOfEdges; ++i) {
        if(!isLifted[i]){
            bool strongBase=pInstance->checkStrongBase(listOfVertices[i],listOfVertices[i+1]);
            isStrongBase[i]=strongBase;
        }
    }
}


double ldp_path_factor::getMinMarginal(const size_t &edgeIndex,const std::vector<double>* pCosts) const{
    //assert(listOfCosts.size()==listOfVertices.size()&&listOfCosts.size()==isLifted.size());
    assert(listOfCosts.size()==numberOfEdges);
    assert(numberOfEdges>=2);
    assert(edgeIndex<numberOfEdges);
    //std::cout<<" printing path factor "<<std::endl;
    for(size_t i=0;i<numberOfEdges;i++){
        //std::cout<<"vertex "<<listOfVertices[i]<<std::endl;
        //std::cout<<"vertex "<<listOfVertices[i]<<std::endl;
    }

    double restrictOne=minimize(pCosts,edgeIndex,1);
    //std::cout<<"restrict one "<<restrictOne<<std::endl;
    double restrictZero=minimize(pCosts,edgeIndex,0);
   // std::cout<<"restrict zero "<<restrictZero<<std::endl;
    return restrictOne-restrictZero;
}


std::vector<double> ldp_path_factor::getAllMinMarginals()const{
    std::vector<double> localCosts=listOfCosts;
    std::vector<double> minMarginals(numberOfEdges);
    double currentOpt=minimize(&localCosts,numberOfEdges,false);
    std::vector<char> origOptSolution=optSolution;
    for (size_t i = 0; i < numberOfEdges; ++i) {
        if(origOptSolution[i]){
            double restrictedOpt=minimize(&localCosts,i,false);
            double delta=currentOpt-restrictedOpt;
            currentOpt=restrictedOpt;
            minMarginals[i]=delta;
            localCosts[i]-=delta;
        }
        else{
            double restrictedOpt=minimize(&localCosts,i,true);
            double delta=restrictedOpt-currentOpt;
            minMarginals[i]=delta;
            localCosts[i]-=delta;
        }
    }
    return minMarginals;
}


void ldp_path_factor::init_primal(){
    std::fill(primalSolution.begin(),primalSolution.end(),0);
}


void ldp_path_factor::print() const{
    std::cout<<"path factor, is must cut "<<mustCut<<std::endl;
    for (size_t i = 0; i < listOfVertices.size(); ++i) {
        char il=isLifted.at(i);
        std::cout<<"vertex "<<listOfVertices[i]<<", is lifted "<<int(il)<<", cost "<<listOfCosts[i]<<std::endl;
    }

}

void ldp_path_factor::updateEdgeCost(const size_t& edgeIndex,const double& value){
    assert(edgeIndex<listOfCosts.size());
   // std::cout<<"update edge cost, index "<<edgeIndex<<", value "<<value<<std::endl;
    listOfCosts[edgeIndex]+=value;
}

void ldp_path_factor::setPrimal(const std::vector<size_t>& primalDescendants, const std::vector<size_t> &vertexLabels){

    size_t numberOfZeros=0;
    size_t indexOfZero=0;
    size_t maxNumber=numberOfEdges;
    if(!mustCut){
        maxNumber--;
    }
    for (size_t i=0;i<listOfVertices.size()-1;i++) {

        size_t vertex1=listOfVertices[i];
        size_t vertex2=listOfVertices[i+1];
        assert(vertex1<vertexLabels.size()&&vertex2<vertexLabels.size()&&vertex1<primalDescendants.size());

        if(isLifted[i]){
            if(vertexLabels[vertex1]==vertexLabels[vertex2]&&vertexLabels[vertex1]!=0){
                primalSolution[i]=1;
            }
            else{
                primalSolution[i]=0;
                numberOfZeros++;
                indexOfZero=i;
            }
        }
        else{
            if(primalDescendants[vertex1]==vertex2){
                primalSolution[i]=1;
            }
            else {
                primalSolution[i]=0;
                indexOfZero=i;
                numberOfZeros++;
            }

        }
    }
    if(!mustCut){
        size_t vertex2=listOfVertices.back();
        size_t vertex1=listOfVertices[0];
        primalSolution.back()=(vertexLabels[vertex1]!=0&&vertexLabels[vertex1]==vertexLabels[vertex2]);
        if(primalSolution.back()==0){
            indexOfZero=numberOfEdges-1;
            numberOfZeros++;
        }
        assert(numberOfZeros!=1||!isLifted.at(indexOfZero));
    }

}

double ldp_path_factor::EvaluatePrimal() const{
    double value=0;
    primalBaseCost=0;
    primalLiftedCost=0;
    for (size_t i=0;i<numberOfEdges;i++) {
        if(primalSolution[i]){
            value+=listOfCosts[i];
            if(isLifted[i]){
                primalLiftedCost+=listOfCosts[i];
            }
            else{
                primalBaseCost+=listOfCosts[i];
            }
        }
    }
    return value;
}


//double ldp_path_factor::minimizeForMustCut(const std::vector<double>*pCosts,size_t edgeIndex,bool edgeLabel)const{
//    const std::vector<double>& costs=*pCosts;
//    double optimalValue=0;
//    size_t indexOfWeakestNegative=numberOfEdges;
//    size_t numberOfPositive=1;
//    size_t indexOfPositive=numberOfEdges;
//    std::fill(optSolution.begin(),optSolution.end(),1);
//    assert(optSolution.size()==numberOfEdges);
//    assert(costs.size()==numberOfEdges);

//    if(edgeIndex<numberOfEdges){
//        assert(edgeIndex<numberOfEdges-1);


//        if(edgeLabel){
//            optimalValue=costs[edgeIndex];
//        }
//        else{
//            numberOfPositive=1;
//            indexOfPositive=edgeIndex;
//            optSolution[indexOfPositive]=0;
//        }
//    }


//    for(size_t i=0;i<numberOfEdges-1;i++){
//        if(i==edgeIndex) continue;
//        if(costs[i]<0){
//            optimalValue+=costs[i];
//            if(indexOfWeakestNegative==numberOfEdges){
//                indexOfWeakestNegative=i;
//            }
//            else if(costs[i]>costs[indexOfWeakestNegative]){
//                indexOfWeakestNegative=i;
//            }
//        }
//        else{
//            numberOfPositive++;
//            optSolution[i]=0;
//            if(numberOfPositive==1){
//                indexOfPositive=i;
//            }
//        }

//    }

//    assert(numberOfPositive==0||indexOfPositive<numberOfEdges);

//    if(numberOfPositive==1&&(isLifted[indexOfPositive]||isStrongBase[indexOfPositive])){
//        assert(indexOfWeakestNegative<numberOfEdges);
//        assert(indexOfPositive<numberOfEdges);
//        if(edgeIndex<numberOfEdges&&!edgeLabel){  //must cut the edge

//            double value=optimalValue-costs[indexOfWeakestNegative];
//            optSolution[indexOfWeakestNegative]=0;
//          //  std::cout<<"restrict zero contradiction, return  "<<value<<std::endl;
//            return value;
//        }
//        else{ //can choose between cut one more edge or join all

//            double allActive=optimalValue+costs[indexOfPositive];
//           // if(edgeIndex<numberOfEdges)std::cout<<"all active cost "<<allActive<<std::endl;
//            double cut=optimalValue-costs[indexOfWeakestNegative];
//           // if(edgeIndex<numberOfEdges)std::cout<<"do cut cost "<<cut<<std::endl;
//            if(cut<allActive){
//                optSolution[indexOfWeakestNegative]=0;
//                return cut;
//            }
//            else{
//                optSolution[indexOfPositive]=1;
//                return allActive;
//            }

//        }
//    }
//    else{  //no contradiction
//       // if(edgeIndex<numberOfEdges) std::cout<<"no cotradiction "<<std::endl;
//        return optimalValue;
//    }


//}





double ldp_path_factor::minimize(const std::vector<double>*pCosts,size_t edgeIndex,bool edgeLabel)const{
    const std::vector<double>& costs=*pCosts;
    double optimalValue=0;
    size_t indexOfWeakestNegative=numberOfEdges;
    size_t numberOfPositive=0;
    if(mustCut){
        numberOfPositive=1;
    }

    size_t indexOfPositive=numberOfEdges;
    std::fill(optSolution.begin(),optSolution.end(),1);
    assert(optSolution.size()==numberOfEdges);
    assert(costs.size()==numberOfEdges);

    if(edgeIndex<numberOfEdges){

        if(edgeLabel){
           // std::cout<<"minimize with restrict one "<<edgeIndex<<std::endl;
            optimalValue=costs[edgeIndex];
        }
        else{
            //std::cout<<"minimize with restrict zero "<<edgeIndex<<std::endl;
            numberOfPositive++;
            indexOfPositive=edgeIndex;
            optSolution[indexOfPositive]=0;
        }
    }


    for(size_t i=0;i<numberOfEdges;i++){
        if(i==edgeIndex) continue;
        if(costs[i]<0){
            optimalValue+=costs[i];
            if(indexOfWeakestNegative==numberOfEdges){
                indexOfWeakestNegative=i;
            }
            else if(costs[i]>costs[indexOfWeakestNegative]){
                indexOfWeakestNegative=i;
            }
        }
        else{
            numberOfPositive++;
            optSolution[i]=0;
            indexOfPositive=i;
        }

    }

    assert(mustCut||numberOfPositive==0||indexOfPositive<numberOfEdges);

    if(numberOfPositive==1&&(mustCut||isLifted[indexOfPositive]||isStrongBase[indexOfPositive])){
        assert(indexOfWeakestNegative<numberOfEdges);
        assert(mustCut||indexOfPositive<numberOfEdges);

        if(mustCut||(edgeIndex<numberOfEdges&&!edgeLabel)){  //must cut the edge

            double value=optimalValue-costs[indexOfWeakestNegative];
            optSolution[indexOfWeakestNegative]=0;
          //  std::cout<<"restrict zero contradiction, return  "<<value<<std::endl;
            return value;
        }
        else{ //can choose between cut one more edge or join all

            double allActive=optimalValue+costs[indexOfPositive];
           // if(edgeIndex<numberOfEdges)std::cout<<"all active cost "<<allActive<<std::endl;
            double cut=optimalValue-costs[indexOfWeakestNegative];
           // if(edgeIndex<numberOfEdges)std::cout<<"do cut cost "<<cut<<std::endl;
            if(cut<allActive){
                optSolution[indexOfWeakestNegative]=0;
                return cut;
            }
            else{
                optSolution[indexOfPositive]=1;
                return allActive;
            }

        }
    }
    else{  //no contradiction
       // if(edgeIndex<numberOfEdges) std::cout<<"no cotradiction "<<std::endl;
        return optimalValue;
    }


}


double ldp_path_factor::LowerBound()const{
    double optValue=minimize(&listOfCosts,numberOfEdges,0);
    if(debug()){
        double value=0;
        for (size_t i = 0; i < numberOfEdges; ++i) {
            if(optSolution[i]) value+=listOfCosts[i];
        }
        assert(std::abs(value-optValue)<eps);
    }
    return optValue;

}

}
