#include "lifted_disjoint_paths/ldp_paths_extractor.hxx"

namespace LPMP {


std::vector<size_t> LdpPathsExtractor::getTimeLayersAfterPaths()const{
    std::vector<size_t> paths(vg.getMaxTime()-maxPathsTime);
    for (size_t i = 0; i < paths.size(); ++i) {
        assert(i+maxPathsTime+1<=vg.getMaxTime());
        paths[i]=vg.getGroupVertices(i+maxPathsTime+1).size();
    }
    return paths;

}

std::vector<size_t> LdpPathsExtractor::getTimeLayersBeforePaths()const{
    std::vector<size_t> paths(minPathsTime-1);
    for (size_t i = 0; i < paths.size(); ++i) {
        paths[i]=vg.getGroupVertices(i+1).size();
    }
    return paths;
}

void LdpPathsExtractor::printExtractedPaths()const{
    for (size_t i = 0; i < extractedPaths.size(); ++i) {
        for (size_t j = 0; j < extractedPaths[i].size(); ++j) {
            std::cout<<extractedPaths[i][j]<<",";
        }
        std::cout<<std::endl<<std::endl;

    }
}


LdpPathsExtractor::LdpPathsExtractor(const VertexGroups<>& vertexGroups,const std::vector<std::vector<size_t>>& paths, size_t cutoff,bool isFirst,bool isLast,size_t _vertexShift):
       vg(vertexGroups)
{
    minIntervalVertex=_vertexShift;
    isFirstInterval=isFirst;
    isLastInterval=isLast;

    if(isFirst){
        minPathsTime=1;
    }
    else{
        minPathsTime=cutoff+1;    }
    if(isLast){
        maxPathsTime=vg.getMaxTime();
    }
    else{
        assert(vg.getMaxTime()>cutoff);
        maxPathsTime=vg.getMaxTime()-cutoff;
    }
    assert(maxPathsTime>minPathsTime);
    extractedPaths=vg.extractInnerPaths(paths,minPathsTime,maxPathsTime,minIntervalVertex);
    minPathsVertex=vg.getMinVertexInTime(minPathsTime)+minIntervalVertex;
    maxPathsVertex=vg.getMaxVertexInTime(maxPathsTime)+minIntervalVertex;
    assert(maxPathsVertex>minPathsVertex);

    vertexToPath=std::vector<size_t>(maxPathsVertex-minPathsVertex+1);
    for (size_t i = 0; i < extractedPaths.size(); ++i) {
        for (size_t j = 0; j < extractedPaths[i].size(); ++j) {
            size_t vertex=extractedPaths[i][j];
            assert(vertex>=minPathsVertex);
            assert(vertex<=maxPathsVertex);
            assert(vertex-minPathsVertex<vertexToPath.size());
            vertexToPath.at(vertex-minPathsVertex)=i;
        }

    }


}

}
