#include "lifted_disjoint_paths/ldp_interval_connection.hxx"


namespace LPMP{

LdpIntervalConnection::LdpIntervalConnection( const LdpPathsExtractor& _pathExtractor1,
                                              const LdpPathsExtractor& _pathExtractor2):
    pathExtractor1(_pathExtractor1),
    pathExtractor2(_pathExtractor2)
{
    constructorBegin=std::chrono::steady_clock::now();
    n1=pathExtractor1.getExtractedPaths().size();
    n3=pathExtractor2.getExtractedPaths().size();
    firstFreeVertex=pathExtractor1.getMaxPathsVertex()+1;
    assert(pathExtractor2.getMinPathsVertex()>1);
    lastFreeVertex=pathExtractor2.getMinPathsVertex()-1;
    assert(firstFreeVertex<lastFreeVertex);
    n2=lastFreeVertex-firstFreeVertex+1;
    numberOfLocalVertices=n1+n2+n3;
    verticesScore=std::vector<double>(numberOfLocalVertices);
    minVertex=pathExtractor1.getMinIntevalVertex();

    std::cout<<"first extracted paths"<<std::endl;
    pathExtractor1.printExtractedPaths();
    std::cout<<"second extracted paths"<<std::endl;
    pathExtractor2.printExtractedPaths();
  //  pResultStructure=nullptr;

}


void LdpIntervalConnection::initScoreOfVertices(const std::vector<size_t> &listOfVertices, const std::vector<double> &costs){
    assert(listOfVertices.size()==costs.size());
    std::fill(verticesScore.begin(),verticesScore.end(),0);
    for (size_t i = 0; i < costs.size(); ++i) {
        size_t vertex=listOfVertices[i];
        if(vertex>=n1&&vertex<n1+n2){
            size_t localVertex=globalToLocalFree(vertex);
            verticesScore[localVertex]=costs[i];
        }
    }
}


void LdpIntervalConnection::initFromFile(const std::string& fileName){
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
            if(v<pathExtractor1.getMinPathsVertex()) continue;
            if(v>pathExtractor2.getMaxPathsVertex()) break;


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

            assert(w>v);
            if(v>pathExtractor2.getMaxPathsVertex()) break;
            if(v<pathExtractor1.getMinPathsVertex()) continue;
            if(w>pathExtractor2.getMaxPathsVertex()) continue;

            edges.push_back({v,w});
            double score = std::stod(strings[2]);
            costs.push_back(score);
        }
       // std::cout<<"edges read"<<std::endl;

        data.close();
       // std::cout<<"init from file done, calling init from vectors"<<std::endl;
        initEdgesFromVectors(edges,costs);
        initScoreOfVertices(vertices,scoreOfVertices);

    }
    catch (std::system_error& er) {
        std::clog << er.what() << " (" << er.code() << ")" << std::endl;

    }

}

void LdpIntervalConnection::initEdgesFromVectors(const std::vector<std::array<size_t, 2> > &edges, const std::vector<double> &costs){
    assert(edges.size()==costs.size());
    std::vector<std::array<size_t,2>> freeEdges;
    std::vector<double> costOfFreeEdges;
    std::vector<std::map<size_t,double>> firstPathsToFree(n1);
    std::map<size_t,std::map<size_t,double>> freeToSecondPaths;
   // std::map<size_t,std::map<size_t,double>> pathsToPaths;

    for (size_t i = 0; i < edges.size(); ++i) {
        size_t vertex1=edges[i][0];
        size_t vertex2=edges[i][1];

        assert(vertex2>vertex1);
        if(vertex1>=pathExtractor1.getMinPathsVertex()&&vertex2<=pathExtractor2.getMaxPathsVertex()){
            char graphPart1=vertexToGraphPartGlobal(vertex1);
            char graphPart2=vertexToGraphPartGlobal(vertex2);
            if(graphPart1==0){
                size_t pathLocalID=firstPathsVertexToLocalIndex(vertex1);
                if(graphPart2==1){
                    size_t vertex2LocalIndex=globalToLocalFree(vertex2);
                    firstPathsToFree[pathLocalID][vertex2LocalIndex]+=costs[i];
                }
                else if(graphPart2==2){
                    if(diagnostics()) std::cout<<"WARNIG: Edges directly between first and second interval tracklets are not supported and are skipped."<<std::endl;
                   // size_t vertex2LocalIndex=secondPathsVertexToLocalIndex(vertex2);
                   // pathsToPaths[pathLocalID][vertex2LocalIndex]+=costs[i];
                }
            }
            else if (graphPart1==1) {
                size_t vertex1LocalIndex=globalToLocalFree(vertex1);
                if(graphPart2==1){
                    size_t vertex2LocalIndex=globalToLocalFree(vertex2);
                    freeEdges.push_back({vertex1LocalIndex,vertex2LocalIndex});
                    costOfFreeEdges.push_back(costs[i]);
                }
                else{
                    assert(graphPart2==2);
                    size_t vertex2LocalIndex=secondPathsVertexToLocalIndex(vertex2);
                    freeToSecondPaths[vertex1LocalIndex][vertex2LocalIndex]+=costs[i];
                }
            }
        }
    }

    for (size_t i = 0; i < n1; ++i) {
        auto edgesFromVertex1=firstPathsToFree[i];
        size_t vertex1=i;
        for (auto iter2=edgesFromVertex1.begin();iter2!=edgesFromVertex1.end();iter2++) {
            size_t vertex2=iter2->first;
            double cost=iter2->second;
            freeEdges.push_back({vertex1,vertex2});
            costOfFreeEdges.push_back(cost);
        }

    }

    for (auto iter=freeToSecondPaths.begin();iter!=freeToSecondPaths.end();iter++) {
        size_t vertex1=iter->first;
        auto edgesFromVertex1=iter->second;
        for (auto iter2=edgesFromVertex1.begin();iter2!=edgesFromVertex1.end();iter2++) {
            size_t vertex2=iter2->first;
            double cost=iter2->second;
            freeEdges.push_back({vertex1,vertex2});
            costOfFreeEdges.push_back(cost);
        }
    }

//    for (auto iter=pathsToPaths.begin();iter!=pathsToPaths.end();iter++) {
//        size_t vertex1=iter->first;
//        auto edgesFromVertex1=iter->second;
//        for (auto iter2=edgesFromVertex1.begin();iter2!=edgesFromVertex1.end();iter2++) {
//            size_t vertex2=iter2->first;
//            double cost=iter2->second;
//            freeEdges.push_back({vertex1,vertex2});
//            costOfFreeEdges.push_back(cost);
//        }
//    }


    EdgeVector edgeVector(freeEdges);
    InfoVector infoVector(costOfFreeEdges);

    completeGraph=LdpDirectedGraph(edgeVector,infoVector,numberOfLocalVertices);

    std::vector<size_t> firstLayers=pathExtractor1.getTimeLayersAfterPaths();
    std::vector<size_t> secondLayers=pathExtractor2.getTimeLayersBeforePaths();

    std::vector<size_t> allLayers;
    allLayers.push_back(n1);
    allLayers.insert(allLayers.end(),firstLayers.begin(),firstLayers.end());
    allLayers.insert(allLayers.end(),secondLayers.begin(),secondLayers.end());
    allLayers.push_back(n3);

    vg.initFromVector(allLayers);

}

char LdpIntervalConnection::vertexToGraphPartGlobal(const size_t& vertexGlobalID)const{
    assert(vertexGlobalID>=pathExtractor1.getMinPathsVertex());
    assert(vertexGlobalID<=pathExtractor2.getMaxPathsVertex());
    if(vertexGlobalID<=pathExtractor1.getMaxPathsVertex()){
        return 0;
    }
    else if(vertexGlobalID<pathExtractor2.getMinPathsVertex()){
        return 1;
    }
    else{
        return 2;
    }

}

char LdpIntervalConnection::vertexToGraphPartLocal(const size_t& localVertexIndex)const{
    assert(localVertexIndex<numberOfLocalVertices);
    if(localVertexIndex<n1){
        return 0;
    }
    else if(localVertexIndex<n1+n2){
        return 1;
    }
    else{
        return 2;
    }
}


void LdpIntervalConnection::createResultsStructures(const std::vector<std::vector<size_t>> &paths){

    middleIntervalPaths=std::vector<std::vector<size_t>>();
    firstToMiddle=std::vector<int> (n1,-1);
    middleToSecond=std::vector<int>();
    // std::vector<size_t> startInMiddle;
    // std::vector<size_t> startInSecond;

    //std::vector<bool>isProcessedInSecond(n3);

    std::vector<std::vector<size_t>> pathsToReturn;
    for (size_t i = 0; i < paths.size(); ++i) {
        std::vector<size_t> newPath;
        size_t firstVertex=paths[i][0];
        char graphPart=vertexToGraphPartLocal(firstVertex);


        if(graphPart!=2&&(graphPart!=0||paths[i].size()>1)){

            //assert(paths[i][1]<n1+n2);
            size_t j=0;
            size_t newPathIndex=middleIntervalPaths.size();
            if(graphPart==0){
                assert(paths[i].size()>1);
                assert(vertexToGraphPartLocal(paths[i][1])!=2);
                firstToMiddle[firstVertex]=int(newPathIndex);
                j=1;
            }

            for (; j < paths[i].size()-1; ++j) {
                size_t vertex=paths[i][j];
                assert(vertexToGraphPartLocal(vertex)==1);
                size_t transformedVertex=localToGlobalFree(vertex);
                assert(transformedVertex>=firstFreeVertex&&transformedVertex<=lastFreeVertex);
                newPath.push_back(transformedVertex);
            }

            char graphPartOfLast=vertexToGraphPartLocal(paths[i].back());
            if(graphPartOfLast==1){
                size_t transformedVertex=localToGlobalFree(paths[i].back());
                newPath.push_back(transformedVertex);
                middleIntervalPaths.push_back(newPath);
                middleToSecond.push_back(-1);

            }else{
                assert(graphPartOfLast==2);
                assert(newPath.size()>0);
                size_t secondPathIndex=secondPathsLocalIndexToPathIndex(paths[i].back());
                middleIntervalPaths.push_back(newPath);
                middleToSecond.push_back(int(secondPathIndex));
            }
        }


    }



}

std::vector<std::vector<size_t>> LdpIntervalConnection::decodePaths(const std::vector<std::vector<size_t>>& paths)const{
    std::vector<std::vector<size_t>> pathsToReturn;
    for (size_t i = 0; i < paths.size(); ++i) {
        std::vector<size_t> newPath;
        size_t firstVertex=paths[i][0];
        char graphPart=vertexToGraphPartLocal(firstVertex);
        if(graphPart==0){
            const std::vector<size_t>& startPath=pathExtractor1.getExtractedPaths().at(firstVertex);
            newPath.insert(newPath.end(),startPath.begin(),startPath.end());
        }
        else if (graphPart==1) {
            newPath.push_back(localToGlobalFree(firstVertex));
        }
        else{
            assert(graphPart==2);
            assert(paths[i].size()==1);
            assert(firstVertex>=n1+n2);
            size_t pathIndex=firstVertex-n1-n2;
            const std::vector<size_t>& startPath=pathExtractor2.getExtractedPaths().at(pathIndex);
            newPath.insert(newPath.end(),startPath.begin(),startPath.end());
        }
        for (size_t j = 1; j < paths[i].size()-1; ++j) {
            size_t vertex=paths[i][j];
            assert(vertexToGraphPartLocal(vertex)==1);
            size_t transformedVertex=localToGlobalFree(vertex);
            assert(transformedVertex>=firstFreeVertex&&transformedVertex<=lastFreeVertex);
            newPath.push_back(transformedVertex);
        }
        if(paths[i].size()>1){
            size_t lastVertex=paths[i].back();
            char graphPart=vertexToGraphPartLocal(lastVertex);
            assert(graphPart!=0);
            if(graphPart==1){
                size_t transformedVertex=localToGlobalFree(lastVertex);
                newPath.push_back(transformedVertex);
            }
            else if(graphPart==2){
                size_t pathIndex=lastVertex-n1-n2;
                const std::vector<size_t>& startPath=pathExtractor2.getExtractedPaths().at(pathIndex);
                newPath.insert(newPath.end(),startPath.begin(),startPath.end());
            }
        }
        pathsToReturn.push_back(newPath);
    }

     for (size_t i = 0; i < pathsToReturn.size(); ++i) {
         for (size_t j = 0; j < pathsToReturn[i].size(); ++j) {
             std::cout<<pathsToReturn[i][j]<<",";
         }
         std::cout<<std::endl<<std::endl;

     }
     return pathsToReturn;
}





}
