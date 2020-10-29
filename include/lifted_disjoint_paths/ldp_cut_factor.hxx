#ifndef LDP_CUT_FACTOR_HXX
#define LDP_CUT_FACTOR_HXX

#endif // LDP_CUT_FACTOR_HXX
#include <utility>
#include<cstdlib>
#include<vector>
#include<array>
#include<set>
#include<map>
#include <config.hxx>

namespace LPMP {


class ldp_cut_factor
{
public:

    ldp_cut_factor(size_t v_,size_t w_, double liftedCost,std::map<size_t,std::map<size_t,double>> inputEdges);
    double LowerBound() const;


    double EvaluatePrimal() const;

    void setPrimal(const std::vector<bool>& primalLabels);

    const std::vector<bool>& getPrimal();

    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }

    auto export_variables() { return std::tie(); }

    void init_primal();

   // std::array<double,3> getAllMinMarginals();

    double getOneMinMarginal(size_t edgeId) const;

    void updateCost(size_t edgeId,double update);

private:
std::vector<double> costs;
std::vector<std::array<size_t,2>> edges;
std::vector<size_t> inputVertices;
std::vector<size_t> outputVertices;
size_t v;
size_t w;
//double liftedEdgeCost;
std::vector<bool> primalSolution;
size_t baseCoveringLifted;





};

ldp_cut_factor::ldp_cut_factor(size_t v_, size_t w_, double liftedCost, std::map<size_t,std::map<size_t,double>> inputEdges):
    v(v_),
    w(w_)
{ //TODO: maybe inputEdges as map<size_t<map<size_t,double>> - better for creating edges as pairs of order of v1 and v2
    std::set<size_t> outVertices;
    size_t edgeCounter=0;
    for(auto iter=inputEdges.begin();iter!=inputEdges.end();iter++){
        size_t v1=iter->first;
        inputVertices.push_back(v1);
        std::map<size_t,double>& neighbors=iter->second;
        for(auto iter2=neighbors.begin();iter2!=neighbors.end();iter2++){
            outVertices.insert(iter2->first);
            edgeCounter++;
        }
    }
    for(size_t vert:outVertices){
        outputVertices.push_back(vert);
    }
    size_t counterInput=0;

    baseCoveringLifted=edgeCounter; //meaning there is none

    for(auto iter=inputEdges.begin();iter!=inputEdges.end();iter++){
        size_t v1=iter->first;
        std::map<size_t,double>& neighbors=iter->second;
        for(size_t i=0;i<outputVertices.size();i++){
            auto it=neighbors.find(outputVertices[i]);
            if(it!=neighbors.end()){
                if(v1==v&&it->first==w){
                    baseCoveringLifted=edges.size();
                }
                edges.push_back({counterInput,i});
                costs.push_back(it->second);
            }
        }
        counterInput++;
    }
    costs.push_back(liftedCost);
    primalSolution=std::vector<bool>(edges.size()+1);
    assert(edgeCounter==edges.size());

}

void ldp_cut_factor::setPrimal(const std::vector<bool>& primalLabels) {
    assert(primalLabels.size()==(edges.size()+1));
    primalSolution=primalLabels;
}

void ldp_cut_factor::updateCost(size_t edgeId,double update){
   costs.at(edgeId)+=update;
}

const std::vector<bool>& ldp_cut_factor::getPrimal(){
    return primalSolution;
}

double ldp_cut_factor::EvaluatePrimal() const{
    double value=0;
    assert(primalSolution.size()==costs.size());
    for (size_t i = 0; i < costs.size(); ++i) {
        value+=primalSolution[i]*costs[i];
    }
    return value;
}

double ldp_cut_factor::LowerBound() const{
    size_t liftedIndex=edges.size();
    double lb=0;
    if(costs.at(liftedIndex)>=0){
        std::vector<bool> usedInput(inputVertices.size(),0);
        std::vector<bool> usedOutput(outputVertices.size(),0);
        bool canUseSimple=true;
        for(size_t i=0;i<costs.size()-1;i++){
            if(costs[i]<0){
                if(i==baseCoveringLifted){
                    if(costs[i]<-costs[liftedIndex]){
                        lb+=costs[i]+costs[liftedIndex];
                        if(usedInput.at(edges[i][0])||usedOutput.at(edges[i][1])){
                            canUseSimple=false;
                            break;
                        }
                        else{
                            usedInput.at(edges[i][0])=1;
                            usedOutput.at(edges[i][1])=1;
                        }
                    }

                }
                else{
                    lb+=costs[i];
                    if(usedInput.at(edges[i][0])||usedOutput.at(edges[i][1])){
                        canUseSimple=false;
                        break;
                    }
                    else{
                        usedInput.at(edges[i][0])=1;
                        usedOutput.at(edges[i][1])=1;
                    }

                }
            }
        }
        if(canUseSimple){
            return lb;
        }
        else{
            //TODO solve min cost assignment
        }
    }
    else{
        std::vector<bool> usedInput(inputVertices.size(),0);
        std::vector<bool> usedOutput(outputVertices.size(),0);
        bool canUseSimple=true;
        double minValue=std::numeric_limits<double>::max();
        for(size_t i=0;i<costs.size()-1;i++){
            minValue=std::min(minValue,costs[i]);
            if(costs[i]<0){
                lb+=costs[i];
                if(usedInput.at(edges[i][0])||usedOutput.at(edges[i][1])){
                    canUseSimple=false;
                    break;
                }
                else{
                    usedInput.at(edges[i][0])=1;
                    usedOutput.at(edges[i][1])=1;
                }
            }
        }
        if(canUseSimple){
            if(minValue<0){
                lb+=costs[liftedIndex];
            }
            else if(costs[liftedIndex]<-minValue){
                lb=costs[liftedIndex]+minValue;
            }
            return lb;
        }
        else{
            //TODO solve min cost assignment
        }

    }

}


double ldp_cut_factor::getOneMinMarginal(size_t edgeId) const{
    assert(edgeId<=costs.size());
    size_t liftedIndex=edges.size();
    std::vector<bool> usedInput(inputVertices.size(),0);
    std::vector<bool> usedOutput(outputVertices.size(),0);
    //std::vector<bool> activeEdges(edges.size(),0); //maybe change to set of indices
    std::set<size_t> activeEdges;
    double lb=0;
    bool canUseSimple=true;
    for(size_t i=0;i<costs.size()-1;i++){  //TODO maybe make one method to be used in LowerBound too!
        if(costs[i]<0){
            lb+=costs[i];
            //activeEdges[i]=1;
            activeEdges.insert(i);
            if(usedInput.at(edges[i][0])||usedOutput.at(edges[i][1])){
                canUseSimple=false;
                break;
            }
            else{
                usedInput.at(edges[i][0])=1;
                usedOutput.at(edges[i][1])=1;
            }
        }
    }
    if(canUseSimple){

        double inactiveLiftedCost=0;
        double activeLiftedCost=0;
        if(costs[baseCoveringLifted]>=0){
            inactiveLiftedCost=lb;
        }
        else{
            inactiveLiftedCost=lb-costs[baseCoveringLifted]; //maybe remove vertices from used input and used output
        }
        if(lb<0){
            activeLiftedCost=lb+costs[liftedIndex];
        }
        else{
            double minValue=costs[0];
            for (size_t i = 1; i < costs.size()-1; ++i) {
                minValue=std::min(minValue,costs[i]);
            }
            activeLiftedCost=minValue+costs[liftedIndex];
        }
        //TODO up to here use for lower bound, return lower from the two
        if(edgeId==liftedIndex){
            return activeLiftedCost-inactiveLiftedCost;
        }
        else if(edgeId!=baseCoveringLifted){
            double optValue=inactiveLiftedCost;
            bool liftedActiveInOpt=false;
            if(activeLiftedCost<inactiveLiftedCost){
                optValue=activeLiftedCost;
                liftedActiveInOpt=true;
            }
            double inactiveCost=0;
            double activeCost=0;
            if(costs[edgeId]>=0){
                inactiveCost=optValue;
                activeCost=costs[edgeId];
                size_t blockV1=edges[edgeId][0];
                size_t blockV2=edges[edgeId][1];
                for(size_t i=0;i<costs.size()-1;i++){
                    if(costs[i]<0&&edges[i][0]!=blockV1&&edges[i][1]!=blockV2){
                        activeCost+=costs[i];
                    }
                }
                if(costs[liftedIndex]<0){
                    activeCost+=costs[liftedIndex];
                }




            }
            else{
                activeCost=optValue;
            }


        }

    }
    else{
        //TODO use min cost assignment
    }

    if(edgeId==liftedIndex){

    }
}


//    ldp_cut_factor(size_t v_,size_t w_, std::vector<std::array<size_t,2>> inputEdges,std::vector<double> inputCosts) { //TODO: maybe inputEdges as map<size_t<map<size_t,double>> - better for creating edges as pairs of order of v1 and v2
//        std::set<size_t> inpVertices;
//        std::set<size_t> outVertices;
//        for (size_t i = 0; i < inputEdges.size(); ++i) {
//            size_t v1=inputEdges[i][0];
//            size_t v2=inputEdges[i][1];
//            if(inpVertices.count(v1)==0){
//                inpVertices.insert(v1);
//            }

//            if(outVertices.count(v2)==0){
//                outVertices.insert(v2);
//            }

//        }
//        //inputVertices=std::vector<size_t>(inpVertices.size());
//        for(size_t vert:inpVertices){
//            inputVertices.push_back(vert);
//        }
//        for(size_t vert:outVertices){
//            outputVertices.push_back(vert);
//        }

//    }



}
