
#include"lifted_disjoint_paths/ldp_cut_factor.hxx"

namespace LPMP{

ldp_cut_factor::ldp_cut_factor(size_t v_, size_t w_, double liftedCost_, std::map<size_t,std::map<size_t,double>> inputEdges):
    v(v_),
    w(w_),
    liftedCost(liftedCost_)
{ //TODO: maybe inputEdges as map<size_t<map<size_t,double>> - better for creating edges as pairs of order of v1 and v2
    std::set<size_t> outVertices;
    numberOfInput=inputEdges.size();
    unassignedLabel=std::numeric_limits<size_t>::max();
    baseCoveringLifted={unassignedLabel,unassignedLabel};
    baseCoverLiftedExists=0;
    liftedActive=false;
    liftedActiveInPrimal=false;
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
    numberOfOutput=outVertices.size();
    size_t counterInput=0;
    std::vector<double> edgeCosts;

    std::vector<std::array<size_t,2>> localEdges(edgeCounter);
    edgeCounter=0;

    for(auto iter=inputEdges.begin();iter!=inputEdges.end();iter++){
        size_t v1=iter->first;
        std::map<size_t,double>& neighbors=iter->second;
        for(size_t i=0;i<outputVertices.size();i++){
            auto it=neighbors.find(outputVertices[i]);
            if(it!=neighbors.end()){
                if(v1==v&&it->first==w){
                    baseCoveringLifted={counterInput,i};
                    baseCoverLiftedExists=1;
                }
                assert(edgeCounter<localEdges.size());
                localEdges[edgeCounter]={counterInput,i};

                edgeCounter++;
                // std::cout<<"adding edge "<<counterInput<<", "<<i<<std::endl;
                edgeCosts.push_back(it->second);
            }
        }
        counterInput++;
    }

    cutGraph=LdpTwoLayerGraph(localEdges,edgeCosts);
    primalSolution=std::vector<size_t>(inputVertices.size()+1,unassignedLabel); //plus one for lifted edge


    storeLabeling=std::vector<size_t>(numberOfInput+1,unassignedLabel);

}

std::pair<LdpTwoLayerGraph, double> ldp_cut_factor::getAllMinMarginals()const{
    double currentOpt=LowerBound();
    bool addLiftedCost=baseCoverLiftedExists&&liftedCost>0;
    LdpTwoLayerGraph localCutGraph=cutGraph;
    LdpTwoLayerGraph minMarginals=cutGraph;

    //double restrictOne= advancedMinimizer(index1,index2,true,addLiftedCost,pCutGraph,pLiftedCost);

    std::vector<size_t> optLabeling=storeLabeling;
    for (size_t i = 0; i < inputVertices.size(); ++i) {
        auto* iter=localCutGraph.forwardNeighborsBegin(i);
        auto *end=localCutGraph.forwardNeighborsEnd(i);
        size_t counter=0;
        for (;iter!=end;iter++) {
            if(counter==optLabeling[i]){
                double restrictZero=advancedMinimizer(i,iter->head,false,addLiftedCost,&localCutGraph,&liftedCost);
                double delta=currentOpt-restrictZero;
                currentOpt=restrictZero;
                iter->cost-=delta;
                minMarginals.setForwardEdgeCost(i,counter,delta);

            }
            else{
                double restrictOne=advancedMinimizer(i,iter->head,true,addLiftedCost,&localCutGraph,&liftedCost);
                double delta=restrictOne-currentOpt;
                iter->cost-=delta;
                minMarginals.setForwardEdgeCost(i,counter,delta);
            }
            counter++;
        }
    }
    double liftedMM=getLiftedMinMarginal(&localCutGraph,&liftedCost);

    return std::pair<LdpTwoLayerGraph,double>(minMarginals,liftedMM);

}



//Returning neighbors indices, not output indices!
//Last entry for lifted edge, different encoding
void ldp_cut_factor::setPrimal(const std::vector<size_t>& primalDescendants, const std::vector<size_t>& vertexLabels) {
    bool setLiftedActive=(vertexLabels[v]!=0&&vertexLabels[v]==vertexLabels[w]);
    for(size_t i=0;i<inputVertices.size();i++){
        const size_t& vertexID=inputVertices[i];
        primalSolution[i]=unassignedLabel;
        const auto* it=cutGraph.forwardNeighborsBegin(i);
        const auto* end=cutGraph.forwardNeighborsEnd(i);
        size_t counter=0;
        for(;it!=end;it++){
            size_t index=it->head;

            assert(index>=0&&index<outputVertices.size());
            size_t neighborID=outputVertices[index];
            if(neighborID==primalDescendants[vertexID]){
                assert(i<primalSolution.size());
                primalSolution[i]=counter;
            }
            counter++;
        }
    }
    assert(v<vertexLabels.size()&&w<vertexLabels.size());

    //Is it a proper encoding?
    if(setLiftedActive){
        liftedActiveInPrimal=true;
        primalSolution.back()=w;
    }
    else{
        liftedActiveInPrimal=false;
        primalSolution.back()=unassignedLabel;
    }
}


double ldp_cut_factor::EvaluatePrimal() const{
    double value=0;
    for(size_t i=0;i<primalSolution.size()-1;i++){
        if(primalSolution[i]!=unassignedLabel){
            value+=cutGraph.getForwardEdgeCost(i,primalSolution[i]);
        }
    }
    primalBaseCost=value;
    primalLiftedCost=0;

    if(primalSolution.back()==w){
        primalLiftedCost=liftedCost;
        value+=liftedCost;  //meaning lifted edge is active
    }
    return value;
}

void ldp_cut_factor::init_primal(){
    for(size_t& v :primalSolution){
        v=unassignedLabel;
    }
}


const std::vector<size_t>& ldp_cut_factor::getPrimal() {
    return primalSolution;
}


void ldp_cut_factor::print()const {
    std::cout<<"CUT "<<v<<", "<<w<<":"<<liftedCost<<std::endl;
    for (size_t i = 0; i < numberOfInput; ++i) {

        auto * it=cutGraph.forwardNeighborsBegin(i);
        auto * end=cutGraph.forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->head;
            double value=it->cost;
            std::cout<<inputVertices[i]<<", "<<outputVertices[outputIndex]<<": "<<value<<std::endl;
        }

    }
}


void ldp_cut_factor::updateCostLifted(const double& value){
    liftedCost+=value;
}



void ldp_cut_factor::updateCostBaseForward(const size_t& inputVertexIndex, const size_t& neighborIndex,const double& value){
    //double oldValue=cutGraph.getForwardEdgeCost(inputVertexIndex,neighborIndex);
    assert(inputVertexIndex<inputVertices.size());
    cutGraph.updateForwardEdgeCost(inputVertexIndex,neighborIndex,value);
}


void ldp_cut_factor::updateCostBaseBackward(const size_t& inputVertexIndex, const size_t& neighborIndex,const double& value){
    //double oldValue=cutGraph.getForwardEdgeCost(inputVertexIndex,neighborIndex);
    assert(inputVertexIndex<outputVertices.size());
    cutGraph.updateBackwardEdgeCost(inputVertexIndex,neighborIndex,value);
}







double ldp_cut_factor::advancedMinimizer(const size_t& index1, const size_t& index2, bool restrictToOne, bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost)const {
    // not restrictToOne .. block only edge
    // restrictToOne ... block both vertices
    //index1 .. index in inputVertices, index2.. index in outputVertices + numberOfInput

    // if(debug()&&index1!=unassignedLabel) std::cout<<"lifted cost "<<(*pLiftedCost)<<", add lifted "<<addLiftedCost<<std::endl;
    LPMP::linear_assignment_problem_input lapInput;
    double minValue=0;
    std::fill(storeLabeling.begin(),storeLabeling.end(),unassignedLabel);
    liftedActive=false;
    bool activeExists=false;
    const double& localLiftedCost=*pLiftedCost;
    const LdpTwoLayerGraph& localCutGraph=*pCutGraph;
    //TODO: Assert of neighborIndex?

    std::tuple<double,size_t,size_t,char> myTuple;
    if(index1!=unassignedLabel&&!restrictToOne){
        // if(debug()) std::cout<<"restric to zero"<<std::endl;
        myTuple=minCutEdge(index1,index2,pCutGraph,pLiftedCost);
    }
    else{
        myTuple=minCutEdge(unassignedLabel,unassignedLabel,pCutGraph,pLiftedCost);
    }
    double minCutValue=std::get<0>(myTuple);
    if(minCutValue>=0){
        // if(debug()&&index1!=unassignedLabel)  std::cout<<"simple method"<<std::endl;
        if(index1!=unassignedLabel&&restrictToOne){ //must join
            // if(debug()&&index1!=unassignedLabel)    std::cout<<"must join"<<std::endl;
            double valueToReturn=0;
            auto* iter=localCutGraph.forwardNeighborsBegin(index1);
            auto* end=localCutGraph.forwardNeighborsEnd(index1);
            bool found=false;
            for(;iter!=end;iter++){
                if(iter->head==index2){
                    valueToReturn=iter->cost;
                    found=true;
                    break;
                }
            }
            assert(found);
            storeLabeling[index1]=index2;
            if(localLiftedCost<0||(index1==baseCoveringLifted[0]&&index2==baseCoveringLifted[1])){
                storeLabeling.back()=w;
                liftedActive=true;
                valueToReturn+=localLiftedCost;
            }
            return valueToReturn;
        }
        else{
            //if(debug()&&index1!=unassignedLabel)    std::cout<<"must cut"<<std::endl;

            double activeCost=localLiftedCost+minCutValue;
            if(activeCost<0){
               // if(debug()&&index1!=unassignedLabel)    std::cout<<"active cost negative"<<std::endl;
                size_t vertex=std::get<1>(myTuple);
                size_t label=std::get<2>(myTuple);
                assert(vertex<storeLabeling.size());
                storeLabeling[vertex]=label;
                storeLabeling.back()=w;
                liftedActive=true;
                return activeCost;
            }
            else{
                //if(debug()&&index1!=unassignedLabel)    std::cout<<"active cost positive"<<std::endl;
                return 0;
            }
        }

    }
    else{
        //if(debug()&&index1!=unassignedLabel)  std::cout<<"ADVANCED METHOD FOR CUT FACTOR"<<std::endl;

        if(index1==unassignedLabel){
            lapInput=createLAStandard(addLiftedCost,pCutGraph,pLiftedCost);
        }
        else if(restrictToOne){
            // if(debug()&&index1!=unassignedLabel)  std::cout<<"restrict to one"<<std::endl;
            assert(index1<inputVertices.size());
            lapInput=laExcludeVertices(index1,index2,addLiftedCost,pCutGraph,pLiftedCost);

            assert(index2<outputVertices.size());
            bool found=false;
            for (auto* it=localCutGraph.forwardNeighborsBegin(index1);it!=localCutGraph.forwardNeighborsEnd(index1);it++) {
                if(it->head==index2){
                    minValue=it->cost;
                    found=true;
                    break;
                }
            }
            assert(found);
            storeLabeling[index1]=index2;
            if(baseCoverLiftedExists&&baseCoveringLifted[0]==index1&&baseCoveringLifted[1]==index2){
                if(addLiftedCost) minValue+=localLiftedCost;
            }
            activeExists=true;

        }
        else{//edge is restricted to zero
            //if(debug()&&index1!=unassignedLabel)  std::cout<<"restrict to zero"<<std::endl;
            assert(index1<inputVertices.size());
            assert(index2<outputVertices.size());
            lapInput=laExcludeEdge(index1,index2,addLiftedCost,pCutGraph,pLiftedCost);
            // std::cout<<"lap input for zero created "<<std::endl;
        }
        MCF::SSP<long,double> mcf(lapInput.no_mcf_nodes(),lapInput.no_mcf_edges());
        lapInput.initialize_mcf(mcf);

        //  std::cout<<"mcf init"<<std::endl;
        minValue+=mcf.solve();


        for (size_t i = 0; i < mcf.no_edges(); ++i) {
            if(mcf.flow(i)>0.99){
                // int label=mcf.head(i)-numberOfInput;
                assert(lapInput.no_left_nodes==numberOfInput);
                size_t label=mcf.head(i)-numberOfInput;
                size_t vertex=mcf.tail(i);
                if(vertex<numberOfInput&&label<numberOfOutput){
                    //  std::cout<<"vertex "<<vertex<<", label "<<label<<std::endl;
                    storeLabeling[vertex]=label;
                    activeExists=true;
                }
            }
        }


        if(baseCoverLiftedExists&&storeLabeling[baseCoveringLifted[0]]==baseCoveringLifted[1]){
            // if(debug()&&index1!=unassignedLabel)  std::cout<<"base cover lifted active"<<std::endl;
            liftedActive=true;
            storeLabeling.back()=w;
            if(!addLiftedCost) minValue+=localLiftedCost;
        }
        else if(activeExists&&localLiftedCost<0){
            // if(debug()&&index1!=unassignedLabel)  std::cout<<"lifted active"<<std::endl;
            minValue+=localLiftedCost;
            liftedActive=true;
            storeLabeling.back()=w;
        }

        // std::cout<<"lifted active "<<liftedActive<<std::endl;

        return minValue;
    }
}



LPMP::linear_assignment_problem_input   ldp_cut_factor::createLAStandard(bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost) const{
    LPMP::linear_assignment_problem_input lapInput;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        auto * it=pCutGraph->forwardNeighborsBegin(i);
        auto * end=pCutGraph->forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->head;
            double value=it->cost;
            if(baseCoverLiftedExists&&addLiftedCost&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                value+=*pLiftedCost;
            }
            if(value<0){
                assert(outputIndex<outputVertices.size());
                // std::cout<<"adding assignment "<<i<<", "<<outputIndex<<": "<<value<<std::endl;
                lapInput.add_assignment(i,outputIndex,value);
            }
        }

    }
    return lapInput;
}


LPMP::linear_assignment_problem_input   ldp_cut_factor::laExcludeEdge(const size_t& v1, const size_t& v2, bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost) const{
    LPMP::linear_assignment_problem_input lapInput;
    bool edgeFound=false;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        auto * it=pCutGraph->forwardNeighborsBegin(i);
        auto * end=pCutGraph->forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->head;
            double value=it->cost;
            if(i==v1&&outputIndex==v2){
                edgeFound=true;
                continue;
            }
            if(baseCoverLiftedExists&&addLiftedCost&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                value+=*pLiftedCost;
            }

            if(value<0){
                assert(outputIndex<outputVertices.size());
                lapInput.add_assignment(i,outputIndex,value);
                // if(debug())std::cout<<"adding assignment "<<i<<", "<<outputIndex<<": "<<value<<std::endl;
            }
        }

    }
    assert(edgeFound);
    return lapInput;
}


LPMP::linear_assignment_problem_input   ldp_cut_factor::laExcludeVertices(const size_t& v1, const size_t& v2, bool addLiftedCost, const LdpTwoLayerGraph *pCutGraph, const double *pLiftedCost) const{
    LPMP::linear_assignment_problem_input lapInput;
    for (size_t i = 0; i < numberOfInput; ++i) {
        lapInput.add_assignment(i,numberOfOutput,0.0);
        if(i==v1){
            continue;
        }
        auto * it=pCutGraph->forwardNeighborsBegin(i);
        auto * end=pCutGraph->forwardNeighborsEnd(i);
        for(;it!=end;it++){
            const size_t& outputIndex=it->head;
            double value=it->cost;
            if(baseCoverLiftedExists&&addLiftedCost&&i==baseCoveringLifted[0]&&outputIndex==baseCoveringLifted[1]){
                value+=*pLiftedCost;
            }
            if(outputIndex==v2) continue;
            if(value<0){
                assert(outputIndex<outputVertices.size());
                lapInput.add_assignment(i,outputIndex,value);
                //  std::cout<<"adding assignment "<<i<<", "<<outputIndex<<": "<<value<<std::endl;
            }
        }


    }
    return lapInput;
}


double ldp_cut_factor::LowerBound() const{
    bool addLifted=(baseCoverLiftedExists&&liftedCost>0);
    double value=advancedMinimizer(unassignedLabel,unassignedLabel,false,addLifted,&cutGraph,&liftedCost);
    if(debug()){
        double controlValue=0;
        for (int i = 0; i < storeLabeling.size()-1; ++i) {
            size_t outputIndex=storeLabeling[i];
            if(outputIndex!=unassignedLabel){
                bool found=false;
                for (auto iter=cutGraph.forwardNeighborsBegin(i);iter!=cutGraph.forwardNeighborsEnd(i);iter++) {
                    if(iter->head==outputIndex){
                        //size_t index2=outputVertices.at(outputIndex);
                        double val=iter->cost;
                        //std::cout<<inputVertices.at(i)<<"->"<<index2<<": "<<val<<std::endl;
                        controlValue+=val;
                        found=true;
                        break;

                    }
                }
                if(!found){
                    throw std::runtime_error("assigned output is not valid");
                }

            }
            //               else{
            //                   std::cout<<inputVertices.at(i)<<" unassigned"<<std::endl;
            //               }

        }
        if(liftedActive){
            controlValue+=liftedCost;
            //std::cout<<"lifted cost "<<liftedCost<<std::endl;
        }
        if(std::abs(controlValue-value)>eps){
            print();
            std::cout<<"control value: "<<controlValue<<", orig value: "<<value<<std::endl;
            throw std::runtime_error("wrong lower bound in cut factor ");
        }
    }


   // std::cout<<"cut factor lb "<<value;
    return value;
}


double ldp_cut_factor::getOneEdgeMinMarginal(const size_t & index1,const size_t & index2,const LdpTwoLayerGraph* pCutGraph,const double* pLiftedCost) const{
    //  if(debug())std::cout<<"base edge min marginal in cut"<<std::endl;
    bool addLiftedCost=baseCoverLiftedExists&&(*pLiftedCost)>0;
    assert(index1<numberOfInput);
    assert(index2<numberOfOutput);
    double restrictOne= advancedMinimizer(index1,index2,true,addLiftedCost,pCutGraph,pLiftedCost);
   // std::cout<<"restrict one "<<restrictOne<<std::endl;
    double restrictZero=advancedMinimizer(index1,index2,false,addLiftedCost,pCutGraph,pLiftedCost);
   // std::cout<<"restrict zero "<<restrictZero<<std::endl;
    return restrictOne-restrictZero;
}


std::tuple<double,size_t,size_t,char> ldp_cut_factor::minCutEdge(size_t index1,size_t index2,const LdpTwoLayerGraph* pCutGraph,const double* pLiftedCost) const{
    double minValue=std::numeric_limits<double>::max();
    size_t v1=0;
    size_t v2=0;
    char cutCoverLiftedNegative=false;
    for(size_t i=0;i<inputVertices.size();i++){
        auto * iter=pCutGraph->forwardNeighborsBegin(i);
        auto * end=pCutGraph->forwardNeighborsEnd(i);
        for(;iter!=end;iter++){
            if(iter->cost<minValue){
                if(index1!=i||index2!=iter->head){
                    minValue=iter->cost;
                    v1=i;
                    v2=iter->head;
                }
            }
            if(baseCoverLiftedExists&&i==baseCoveringLifted[0]&&iter->head==baseCoveringLifted[1]){
                cutCoverLiftedNegative=(iter->cost<0);
            }
        }

    }

    std::tuple<double,size_t,size_t,char>t (minValue,v1,v2,cutCoverLiftedNegative);
    return t;
}


double ldp_cut_factor::getLiftedMinMarginal(const LdpTwoLayerGraph* pCutGraph,const double* pLiftedCost) const{
  //  if(debug()) std::cout<<"lifted min marginal"<<std::endl;
    double localLiftedCost=*pLiftedCost;
    std::tuple<double,size_t,size_t,char>t=minCutEdge(unassignedLabel,unassignedLabel,pCutGraph,pLiftedCost);
    bool baseCoverNegative=std::get<3>(t)>0;
    double minValue=std::get<0>(t);

    if(!baseCoverLiftedExists||!baseCoverNegative){

        if(minValue>=0){
            //if(debug()) std::cout<<"no base cover, minValue positive"<<std::endl;
            double value=localLiftedCost+minValue;
            return value;
        }
        else{
            //  if(debug()) std::cout<<"no base cover, minValue negative"<<std::endl;
            return localLiftedCost;
        }
    }
    else{  //Base cover lifted exists and is negative
        assert(minValue<0);
        //   if(debug()) std::cout<<"base cover negative"<<std::endl;
        double zeroLiftedCost=0.0;
        double lowerBoundNoLift=advancedMinimizer(unassignedLabel,unassignedLabel,false,false,pCutGraph,&zeroLiftedCost); //some cut edge must be active, not effective! Min cut edge searched again!
        double restrictedOne=0;
        double restrictedZero=0;
        assert(lowerBoundNoLift<0);
        if(storeLabeling[baseCoveringLifted[0]]==baseCoveringLifted[1]){
            //if(debug()) std::cout<<"base cover active, lower bound "<<lowerBoundNoLift<<std::endl;

            restrictedOne=lowerBoundNoLift+localLiftedCost;
            restrictedZero=advancedMinimizer(baseCoveringLifted[0],baseCoveringLifted[1],false,false,pCutGraph,&zeroLiftedCost);

            //if(debug()) std::cout<<"restrict zero value "<<restrictedZero<<std::endl;
            return restrictedOne-restrictedZero;
        }
        else{
            //  if(debug())std::cout<<"base cover not active"<<std::endl;
          //  if(debug()) std::cout<<"lb negative"<<std::endl;
            return localLiftedCost;
        }
    }
}

}
