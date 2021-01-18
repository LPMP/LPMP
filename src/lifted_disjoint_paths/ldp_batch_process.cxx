#include "lifted_disjoint_paths/ldp_batch_process.hxx"

namespace LPMP {

LdpBatchProcess::LdpBatchProcess(VertexGroups<>& shiftedGroups, std::vector<std::array<size_t,2>>& shiftedVertexLabels_, size_t maxLabelSoFar_, size_t maxTimeForLabeled_, size_t minTimeToUse, size_t maxTimeToUse){ //shifted vertex labels assume that vertex with zero index is minVertex in VG
  constructorBegin=std::chrono::steady_clock::now();
    size_t unset=std::numeric_limits<size_t>::max();
    pvg=&shiftedGroups;
    if(minTimeToUse!=unset){
        assert(maxTimeToUse!=unset);
        minTime=minTimeToUse;
        maxTime=maxTimeToUse;
        minVertex=pvg->getMinVertexInTime(minTime);
        maxVertex=pvg->getMaxVertexInTime(maxTime);

    }
    else{
        maxVertex=pvg->getMaxVertex();
        minVertex=pvg->getMinVertex();
        minTime=pvg->getMinTime();
        maxTime=pvg->getMaxTime();

    }
   // std::cout<<"min max time and vertices set"<<std::endl;
    numberOfVerticesInBatch=maxVertex-minVertex+1;
    indexToDel=shiftedVertexLabels_.size();

   // origVerticesScore=std::vector<double>(numberOfVerticesInBatch);
    maxTimeForLabeled=maxTimeForLabeled_;
    assert(minTime<=maxTimeForLabeled+1); //Allow max time for labeled zero in case of first layer and no labels
    assert(maxTime>maxTimeForLabeled); //At least one new layer allowed
    minValidVertex=pvg->getMinVertexInTime(maxTimeForLabeled+1);
    if(minValidVertex>minVertex){
        //std::cout<<"initializing labels"<<std::endl;
        shiftedLabels=std::vector<size_t>(minValidVertex-minVertex);
        for (size_t i = 0; i < shiftedVertexLabels_.size(); ++i) {
            size_t vertex=shiftedVertexLabels_[i][0];
            if(vertex<minVertex) continue;
            if(vertex>=minValidVertex){
                indexToDel=i;
                break;
            }
            shiftedLabels[vertex-minVertex]=shiftedVertexLabels_[i][1];
        }

    }

    maxLabelSoFar=maxLabelSoFar_;

   // assert(shiftedLabels.size()>=minValidVertex);
    numberOfUsedLabels=0;
    numberOfOutputVertices=0;
    edgesCreated=false;
    labelsDecoded=false;
  //  std::cout<<"constructor finished "<<std::endl;


}

void LdpBatchProcess::createLocalVG(LPMP::VertexGroups<>& localVG){
    {
            std::unordered_map<size_t,std::vector<size_t>> groups;
            std::vector<size_t> vToGroup;
            size_t layerCounter=1;
             size_t lastVertex=0;
            if(numberOfUsedLabels>0){
                groups[1]=std::vector<size_t>(numberOfUsedLabels);
                for (size_t i = 0; i < numberOfUsedLabels; ++i) {
                    groups[1][i]=i;
                    vToGroup.push_back(1);
                    lastVertex=i;
                }
                layerCounter=2;
            }




            std::vector<size_t> currenGroup;
            for (size_t i = maxTimeForLabeled+1; i <= maxTime; ++i) {
                const std::vector<size_t>& verticesGlobal= pvg->getGroupVertices(i);
                for (size_t j = 0; j < verticesGlobal.size(); ++j) {
                    size_t globalID=verticesGlobal.at(j);
                    size_t localID=globalIndexToLocalIndex(globalID);
                    assert(localID>=numberOfUsedLabels);
                    currenGroup.push_back(localID);
                    lastVertex=localID;
                    assert(vToGroup.size()==localID);
                    vToGroup.push_back(layerCounter);
                }
                groups[layerCounter]=currenGroup;
                currenGroup=std::vector<size_t>();
                layerCounter++;
            }
            size_t s=lastVertex+1;
            size_t t=lastVertex+2;
            std::vector<size_t> sGroup={s};
            std::vector<size_t> tGroup={t};
            groups[layerCounter]=tGroup;
            groups[0]=sGroup;
            vToGroup.push_back(0);
            vToGroup.push_back(layerCounter);
            localVG=VertexGroups(groups,vToGroup);

        }

}



LdpBatchProcess::LdpBatchProcess(VertexGroups<>& shiftedGroups, std::vector<std::array<size_t,2>>& shiftedVertexLabels_, size_t maxLabelSoFar_, size_t maxTimeForLabeled_){ //shifted vertex labels assume that vertex with zero index is minVertex in VG
    size_t unset=std::numeric_limits<size_t>::max();
    LdpBatchProcess(shiftedGroups,shiftedVertexLabels_,maxTimeForLabeled_,maxLabelSoFar_,unset,unset);
}


void LdpBatchProcess::decode(const std::vector<std::vector<size_t>>& paths){
    size_t globalLabelCounter=maxLabelSoFar+1;
    decodedLabels=std::vector<std::array<size_t,2>>();
    std::vector<size_t> localLabels(numberOfOutputVertices);
    for (size_t i = 0; i < paths.size(); ++i) {
        size_t firstVertex=paths.at(i).at(0);
        size_t labelToUse=0;
        if(firstVertex<numberOfUsedLabels){
            labelToUse=localIndexToGlobalLabel[firstVertex];
        }
        else{
            labelToUse=globalLabelCounter;
            globalLabelCounter++;
        }
        for (size_t j = 0; j < paths.at(i).size(); ++j) {
            localLabels[paths.at(i).at(j)]=labelToUse;
        }
    }

    for (int i = numberOfUsedLabels; i < numberOfOutputVertices; ++i) {
        size_t globalID=localIndexToGlobalIndex(i);
        size_t label=localLabels[i];
        if(label>0){
            decodedLabels.push_back({globalID,label});
        }
    }

    maxLabelSoFar=globalLabelCounter-1;
    labelsDecoded=true;


}

void LdpBatchProcess::initVertexScoreFromVector(const std::vector<size_t>& vertexList,const std::vector<double>& costs){
    assert(vertexList.size()==costs.size());
    assert(edgesCreated);
    assert(outputVerticesScore.size()==numberOfOutputVertices);
    size_t i=0;
    while(i<vertexList.size()){
       // std::cout<<"processing vertex "<<i<<std::endl;
        size_t vertex=vertexList[i];
        if(vertex<minValidVertex){
            i++;
            continue;
        }
        if(vertex>maxVertex) break;
        size_t localIndex=globalIndexToLocalIndex(vertex);
        assert(localIndex<numberOfOutputVertices);
        outputVerticesScore[localIndex]=costs[i];

    }
    //std::cout<<"end of processing vertices"<<std::endl;


}

size_t LdpBatchProcess::globalIndexToLocalIndex(const size_t& globalIndex){  //only for vertices from valid layers
    assert(globalIndex>=minValidVertex&&globalIndex<=maxVertex);
    //assert(numberOfUsedLabels<std::numeric_limits<size_t>::max());
    size_t toReturn=globalIndex+numberOfUsedLabels-minValidVertex;
    //std::cout<<"to return "<<toReturn<<std::endl;
   // assert(localIndexToGlobalIndex(toReturn)==globalIndex);
    return toReturn;
}

size_t LdpBatchProcess::localIndexToGlobalIndex(const size_t &localIndex){
    assert(localIndex>=numberOfUsedLabels&&localIndex<numberOfOutputVertices);
    size_t toReturn=localIndex+minValidVertex-numberOfUsedLabels;
   // assert(globalIndexToLocalIndex(toReturn)==localIndex);
    return toReturn;
}

void LdpBatchProcess::initEdgesFromVector(const std::vector<std::array<size_t,2>>& edges,const std::vector<double> &costs){

    size_t i = 0;
    std::vector<std::array<size_t,2>> outputEdges;

    assert(edges.size()==costs.size());
    while (i < edges.size()) {
        size_t v0=edges[i][0];
        size_t v1=edges[i][1];
        assert(v0<v1);
        if(v0<minVertex){
            i++;
            continue;
        }
        if(v0>=minValidVertex) break;
        if(v1>maxVertex){
            i++;
            continue;
        }
        size_t time0=pvg->getGroupIndex(v0);
        size_t time1=pvg->getGroupIndex(v1);
        assert(time0<time1);
        if(time0>maxTimeForLabeled) break;
        if(time1<=maxTimeForLabeled){
            i++;
            continue;
        }

        size_t vertexIndexInLabels=v0-minVertex;
        assert(vertexIndexInLabels<shiftedLabels.size());
        size_t label=shiftedLabels[vertexIndexInLabels];
        if(label>0){
            edgesFromLabeled[label][v1]+=costs.at(i);
        }
        i++;

    }
    numberOfUsedLabels=edgesFromLabeled.size();
   // std::cout<<"map finished, used labels "<<numberOfUsedLabels<<std::endl;
    localIndexToGlobalLabel=std::vector<size_t>(numberOfUsedLabels);
    numberOfOutputVertices=numberOfUsedLabels+maxVertex-minValidVertex+1;
   // std::cout<<"number of output vertices "<<numberOfOutputVertices<<std::endl;
    outputGraph=andres::graph::Digraph<>(numberOfOutputVertices);
    outputVerticesScore=std::vector<double>(numberOfOutputVertices);


    size_t lCounter=0;
    for (auto iter=edgesFromLabeled.begin();iter!=edgesFromLabeled.end();iter++) {
       // std::cout<<"in for"<<std::endl;
        size_t label=iter->first;
        auto neighbors=iter->second;
        localIndexToGlobalLabel.at(lCounter)=label;
        for(auto iter2=neighbors.begin();iter2!=neighbors.end();iter2++){
            size_t vertexGlobID=iter2->first;
            double cost=iter2->second;
            size_t vertexLocalID=globalIndexToLocalIndex(vertexGlobID);
            outputGraph.insertEdge(lCounter,vertexLocalID);
            outputEdges.push_back({lCounter,vertexLocalID});
            outputEdgeCosts.push_back(cost);
        }


        lCounter++;
    }
    assert(lCounter==numberOfUsedLabels);
    for(;i<edges.size();i++){
        size_t v0=edges[i][0];
        size_t v1=edges[i][1];
    //    std::cout<<v0<<","<<v1<<std::endl;
        assert(v0<v1);
        assert(v0>=minValidVertex);
     //   std::cout<<"after assert"<<std::endl;

        if(v0>maxVertex) break;
        if(v1>maxVertex) continue;

        size_t transformedV0=globalIndexToLocalIndex(v0);
        size_t transformedV1=globalIndexToLocalIndex(v1);
        assert(transformedV0<numberOfOutputVertices);
//        if(transformedV1>=numberOfOutputVertices){
//            std::cout<<"tr v1 "<<transformedV1<<", graph size "<<numberOfOutputVertices<<std::endl;
//        }
        assert(transformedV1<numberOfOutputVertices);
      //  if(transformedV1>=1356) std::cout<<"transformed "<<transformedV0<<", "<<transformedV1<<std::endl;
        outputGraph.insertEdge(transformedV0,transformedV1);
        outputEdges.push_back({transformedV0,transformedV1});
        outputEdgeCosts.push_back(costs.at(i));

    }
    outputVerticesScore=std::vector<double>(numberOfOutputVertices);
    EdgeVector ev(outputEdges);
    InfoVector iv(outputEdgeCosts);

    myOutputGraph=LdpDirectedGraph(ev,iv);
    //std::cout<<"edges finished "<<outputGraph.numberOfEdges()<<std::endl;
    edgesCreated=true;

}

void LdpBatchProcess::initFromFile(std::string fileName){

    char delim=',';
    std::string line;
    std::ifstream data;
    try{
        data.open(fileName);
        if(!data){
            throw std::system_error(errno, std::system_category(), "failed to open graph file "+fileName);
        }

        std::getline(data, line);

        std::vector<std::string> strings;
        std::vector<size_t> vertices;
        std::vector<double> scoreOfVertices;

        //Vertices that are not found have score=0. Appearance and disappearance cost are read here.
        while (std::getline(data, line) && !line.empty()) {
            strings = split(line, delim);

            unsigned int v = std::stoul(strings[0]);
            if(v<minVertex) continue;
            if(v>maxVertex) break;


            double c = std::stod(strings[1]);
            vertices.push_back(v);
            scoreOfVertices.push_back(c);

        }
      //  std::cout<<"vertices read"<<std::endl;

        std::vector<std::array<size_t,2>> edges;
        std::vector<double> costs;

        while (std::getline(data, line) && !line.empty()) {


            strings = split(line, delim);

            unsigned int v = std::stoul(strings[0]);

            unsigned int w = std::stoul(strings[1]);

            if(v>maxVertex) break;
            if(v<minVertex) continue;
            if(w>maxVertex) continue;

            edges.push_back({v,w});
            double score = std::stod(strings[2]);
            costs.push_back(score);
        }
       // std::cout<<"edges read"<<std::endl;

        data.close();
       // std::cout<<"init from file done, calling init from vectors"<<std::endl;
        initEdgesFromVector(edges,costs);
        initVertexScoreFromVector(vertices,scoreOfVertices);

    }
    catch (std::system_error& er) {
        std::clog << er.what() << " (" << er.code() << ")" << std::endl;

    }

}


}
