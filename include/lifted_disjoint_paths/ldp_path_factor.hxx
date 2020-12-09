#ifndef LDP_PATH_FACTOR_HXX
#define LDP_PATH_FACTOR_HXX

#include<cstdlib>
#include<vector>
#include <config.hxx>
#include"lifted_disjoint_paths/ldp_instance.hxx"

namespace LPMP {

class ldp_path_factor{

public:
    ldp_path_factor(const std::vector<size_t>& _listOfVertices,const std::vector<double>& _listOfCosts,const std::vector<char>& _isLifted,const lifted_disjoint_paths::LdpInstance* pInstance):
        listOfVertices(_listOfVertices),
        listOfCosts(_listOfCosts),
        isLifted(_isLifted)
    {
        assert(listOfCosts.size()==listOfVertices.size()&&listOfCosts.size()==isLifted.size());
        numberOfEdges=listOfCosts.size();
        isStrongBase=std::vector<char>(numberOfEdges);
        primalSolution=std::vector<char>(numberOfEdges);
        optSolution=std::vector<char>(numberOfEdges);
        assert(listOfCosts.size()>2);
        for (size_t i = 0; i < numberOfEdges-1; ++i) {
            if(!isLifted[i]){
                bool strongBase=pInstance->isStrongBase(listOfVertices[i],listOfVertices[i+1]);
                isStrongBase[i]=strongBase;
            }
        }
    }



    double LowerBound() const;

    double EvaluatePrimal() const;

    void setPrimal(const std::vector<size_t>& primalDescendants, const std::vector<size_t> &vertexLabels);

    const std::vector<char>& getPrimal(){
        return primalSolution;
    }

    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(listOfCosts);}

    auto export_variables() { return std::tie(listOfCosts); }

    void init_primal();

    // std::array<double,3> getAllMinMarginals();


    void updateEdgeCost(const size_t& edgeIndex,const double& value);

    double getMinMarginal(const size_t& edgeIndex,const std::vector<double> *pCosts) const;


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


    const double& getPrimalLiftedCost()const{
        return primalLiftedCost;
    }
    const double& getPrimalBaseCost()const{
        return primalBaseCost;
    }

    std::vector<double> getAllMinMarginals()const;

    const void print() const;

private:
    double minimize(const std::vector<double> *pCosts, size_t edgeIndex, bool edgeLabel)const;

    const std::vector<size_t> listOfVertices;
    std::vector<double> listOfCosts;
    const std::vector<char> isLifted;
    std::vector<char> isStrongBase;
    size_t numberOfEdges;
    std::vector<char> primalSolution;
    mutable double primalBaseCost;
    mutable double primalLiftedCost;
    mutable std::vector<char> optSolution;

};


double ldp_path_factor::getMinMarginal(const size_t &edgeIndex,const std::vector<double>* pCosts) const{
    assert(listOfCosts.size()==listOfVertices.size()&&listOfCosts.size()==isLifted.size());
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
    for (int i = 0; i < numberOfEdges; ++i) {
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


const void ldp_path_factor::print() const{
    std::cout<<"path factor"<<std::endl;
    for (int i = 0; i < listOfVertices.size(); ++i) {
        int il=isLifted.at(i);
        std::cout<<"vertex "<<listOfVertices[i]<<", is lifted "<<il<<", cost "<<listOfCosts[i]<<std::endl;
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
    for (size_t i=0;i<numberOfEdges-1;i++) {

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
    size_t vertex2=listOfVertices.back();
    size_t vertex1=listOfVertices[0];
    primalSolution.back()=(vertexLabels[vertex1]!=0&&vertexLabels[vertex1]==vertexLabels[vertex2]);
    if(primalSolution.back()==0){
        indexOfZero=numberOfEdges-1;
        numberOfZeros++;
    }
     assert(numberOfZeros!=1||!isLifted.at(indexOfZero));

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

double ldp_path_factor::minimize(const std::vector<double>*pCosts,size_t edgeIndex,bool edgeLabel)const{
    const std::vector<double>& costs=*pCosts;
    double optimalValue=0;
    size_t indexOfWeakestNegative=numberOfEdges;
    size_t numberOfPositive=0;
    size_t indexOfPositive=numberOfEdges;
    std::fill(optSolution.begin(),optSolution.end(),1);
    assert(optSolution.size()==numberOfEdges);
    assert(costs.size()==numberOfEdges);
    //std::cout<<"calling minimize, printing edge costs "<<std::endl;
//    for(size_t i=0;i<numberOfEdges;i++){
//        double c=costs[i];
//        size_t vertex=listOfVertices[i];
//       // std::cout<<listOfVertices[i]<<": "<<c<<std::endl;

//    }
    if(edgeIndex<numberOfEdges){
        //        std::cout<<"printing edge costs "<<std::endl;
        //        for(auto& c:listOfCosts){
        //            std::cout<<c<<std::endl;
        //        }

        if(edgeLabel){
           // std::cout<<"minimize with restrict one "<<edgeIndex<<std::endl;
            optimalValue=costs[edgeIndex];
        }
        else{
            //std::cout<<"minimize with restrict zero "<<edgeIndex<<std::endl;
            numberOfPositive=1;
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
            if(numberOfPositive==1){
                indexOfPositive=i;
            }
        }

    }

    assert(numberOfPositive==0||indexOfPositive<numberOfEdges);

    if(numberOfPositive==1&&(isLifted[indexOfPositive]||isStrongBase[indexOfPositive])){
        assert(indexOfWeakestNegative<numberOfEdges);
        assert(indexOfPositive<numberOfEdges);
        if(edgeIndex<numberOfEdges&&!edgeLabel){  //must cut the edge

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


class ldp_snc_path_message
{
public:
    ldp_snc_path_message(const std::vector<size_t>& _edgeIndexInPath,  //Mostly one, two for the first and the last path vertices
                         const std::vector<size_t>& _vertexIndexInSnc,
                         const std::vector<char>& _isLifted, bool debugInfo=false)
    {
        edgeIndexInPath=_edgeIndexInPath;  //Mostly one, two for the first and the last path vertices
        vertexIndexInSnc=_vertexIndexInSnc;
        isLifted=_isLifted;
        dimension=_edgeIndexInPath.size();
        debInfo=debugInfo;

        assert(dimension==1||dimension==2);
        assert(vertexIndexInSnc.size()==dimension);
        assert(isLifted.size()==dimension);
        assert(edgeIndexInPath.size()==dimension);

    }

    template<typename PATH_FACTOR>
    void RepamLeft(PATH_FACTOR& l, const double msg, const std::size_t msg_dim) const
    {

        assert(msg_dim<dimension);
        assert(vertexIndexInSnc.size()==dimension);
        assert(isLifted.size()==dimension);
        assert(edgeIndexInPath.size()==dimension);

        if(debug()){
            size_t indexInPath=edgeIndexInPath[msg_dim];
            size_t v0=l.getListOfVertices().at(indexInPath);
            size_t v1;
            if(indexInPath==l.getNumberOfEdges()-1){
                v1=l.getListOfVertices().at(0);
                lastV1=v1;
                lastV2=v0;
                lastValue=msg;
            }
            else{
                v1=l.getListOfVertices().at(indexInPath+1);
                lastV1=v0;
                lastV2=v1;
                lastValue=msg;
            }
            assert(isLifted.at(msg_dim)==l.getLiftedInfo().at(indexInPath));
            //std::cout<<"Update cost of path, vertices "<<v0<<", "<<v1<<", edge index "<<edgeIndexInPath[msg_dim]<<", value "<<msg<<std::endl;
        }

        l.updateEdgeCost(edgeIndexInPath[msg_dim],msg);


    }

    template<typename SINGLE_NODE_CUT_FACTOR>
    void RepamRight(SINGLE_NODE_CUT_FACTOR& r, const double msg, const std::size_t msg_dim) const
    {

        assert(dimension==1||dimension==2);
        assert(msg_dim<dimension);
        assert(vertexIndexInSnc.size()==dimension);
        assert(isLifted.size()==dimension);
        assert(edgeIndexInPath.size()==dimension);

        r.updateEdgeCost(msg,vertexIndexInSnc[msg_dim],isLifted[msg_dim]);
        if(debug()){
            size_t centralNodeID=r.nodeID;

            size_t secondVertex;
            size_t v0;
            size_t v1;

            if(isLifted.at(msg_dim)){

                secondVertex=r.getLiftedIDs().at(vertexIndexInSnc[msg_dim]);

            }
            else{
                secondVertex=r.getBaseIDs().at(vertexIndexInSnc[msg_dim]);
            }
            v0=std::min(centralNodeID,secondVertex);
            v1=std::max(centralNodeID,secondVertex);
            assert(v0==lastV1);
            assert(v1==lastV2);
            if(std::abs(msg+lastValue)>=eps){
                std::cout<<"WRONG update cost "<<v0<<", "<<v1<<", left was "<<lastValue<<", right is "<<msg<<std::endl;
            }
            assert(std::abs(msg+lastValue)<eps);
           // std::cout<<"Update cost of SNC, central node "<<centralNodeID<<", second vertex "<< secondVertex<<", value "<<msg<<std::endl;
        }


    }

    template<typename SINGLE_NODE_CUT_FACTOR, typename MSG>
    void send_message_to_left(const SINGLE_NODE_CUT_FACTOR& r, MSG& msg, const double omega = 1.0)
    {
        assert(dimension==1||dimension==2);
        assert(vertexIndexInSnc.size()==dimension);
        assert(isLifted.size()==dimension);
        assert(edgeIndexInPath.size()==dimension);
        assert(isLifted.size()==dimension);

        double delta=0;
        double controlDelta=0;
        if(dimension==1){
            const std::vector<double>& baseCosts=r.getBaseCosts();
            const std::vector<double>& liftedCosts=r.getLiftedCosts();
            if(isLifted.at(0)){
                delta = r.getOneLiftedMinMarginal(vertexIndexInSnc[0],&baseCosts,&liftedCosts);

                if(debug()){
                    std::vector<double> controlLiftedCosts=liftedCosts;
                    controlLiftedCosts[vertexIndexInSnc[0]]-=delta;
                    controlDelta=r.getOneLiftedMinMarginal(vertexIndexInSnc[0],&baseCosts,&controlLiftedCosts);
                    if(abs(controlDelta)>eps){
                        throw std::runtime_error("wrong control snc lifted min marginal in path message");
                    }

                }
            }
            else{

                delta = r.getOneBaseEdgeMinMarginal(vertexIndexInSnc[0],&baseCosts,&liftedCosts);

                if(debug()){
                    std::vector<double> controlBaseCosts=baseCosts;
                    controlBaseCosts[vertexIndexInSnc[0]]-=delta;
                    controlDelta=r.getOneBaseEdgeMinMarginal(vertexIndexInSnc[0],&controlBaseCosts,&liftedCosts);
                    if(abs(controlDelta)>eps){
                        throw std::runtime_error("wrong control snc base min marginal in path message");
                    }
                }
            }

            msg[0] -= omega * delta;
        }
        else{
            std::vector<double> baseCosts=r.getBaseCosts();
            std::vector<double> liftedCosts=r.getLiftedCosts();
            for (size_t i=0;i<dimension;i++) {
                if(isLifted.at(i)){
                   // std::cout<<"lifted "<<std::endl;
                    delta = r.getOneLiftedMinMarginal(vertexIndexInSnc[i],&baseCosts,&liftedCosts);
                    liftedCosts[vertexIndexInSnc[i]]-=delta;

                    if(debug()){
                        controlDelta=r.getOneLiftedMinMarginal(vertexIndexInSnc[i],&baseCosts,&liftedCosts);
                        if(abs(controlDelta)>eps){
                            throw std::runtime_error("control snc min marginal unequall to zero");
                        }
                    }

                }
                else{
                    //std::cout<<"base "<<std::endl;
                    delta = r.getOneBaseEdgeMinMarginal(vertexIndexInSnc[i],&baseCosts,&liftedCosts);
                    baseCosts[vertexIndexInSnc[i]]-=delta;
                    if(debug()){
                        controlDelta=r.getOneBaseEdgeMinMarginal(vertexIndexInSnc[i],&baseCosts,&liftedCosts);
                        if(abs(controlDelta)>eps){
                            throw std::runtime_error("control snc min marginal unequall to zero");
                        }
                    }

                }
                msg[i] -= omega * delta;
            }
        }

    }

    template<typename PATH_FACTOR, typename MSG>
    void send_message_to_right(const PATH_FACTOR& l, MSG& msg, const double omega)
    {

        double delta=0;
        double controlDelta=0;
        if(dimension==1){
            const std::vector<double>& costs=l.getCosts();
            delta=l.getMinMarginal(edgeIndexInPath[0],&costs);
            if(debug()){
                std::vector<double> controlCosts=costs;
                controlCosts.at(edgeIndexInPath[0])-=delta;
                controlDelta=l.getMinMarginal(edgeIndexInPath[0],&controlCosts);
                if(abs(controlDelta)>eps){
                    std::cout<<"control value from path factor"<<controlDelta<<", original delta: "<<delta<<std::endl;
                    throw std::runtime_error("wrong path min marginal check in path message");
                }
            }
            msg[0] -= omega * delta;
        }
        else{
            std::vector<double> localCosts=l.getCosts();

            for (size_t i=0;i<dimension;i++) {
                delta=l.getMinMarginal(edgeIndexInPath[i],&localCosts);
              //  std::cout<<"sending "<<delta<<std::endl;
                assert(edgeIndexInPath[i]<localCosts.size());
                localCosts[edgeIndexInPath[i]]-=delta;
                if(debug()){
                    controlDelta=l.getMinMarginal(edgeIndexInPath[i],&localCosts);
                    if(abs(controlDelta)>eps){
                        std::cout<<"control value from path factor"<<controlDelta<<", original delta: "<<delta<<std::endl;
                        throw std::runtime_error("wrong path min marginal check in path message");
                    }
                }
                msg[i] -= omega * delta;
            }
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
    bool debInfo;

    mutable size_t lastV1;
    mutable size_t lastV2;
    mutable double lastValue;

};




}


#endif // LDP_PATH_FACTOR_HXX
