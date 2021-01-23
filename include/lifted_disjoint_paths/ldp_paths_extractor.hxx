#ifndef LDP_PATHS_EXTRACTOR_HXX
#define LDP_PATHS_EXTRACTOR_HXX

#include "ldp_vertex_groups.hxx"
#include "config.hxx"

namespace LPMP {

class LdpPathsExtractor{
public:
    LdpPathsExtractor( const VertexGroups<>& vertexGroups,const std::vector<std::vector<size_t>>& paths, size_t cutoff,bool isFirst,bool isLast,size_t _vertexShift);
//    const VertexGroups<>& getVertexGroups()const{
//        return vg;
//    }


    const size_t & getMinPathsVertex()const{
        return minPathsVertex;
    }

    const size_t & getMaxPathsVertex()const{
        return maxPathsVertex;
    }

    const size_t & getMinIntevalVertex()const{
        return minIntervalVertex;
    }


    const std::vector<std::vector<size_t>>& getExtractedPaths()const{
        return extractedPaths;
    }

    void printExtractedPaths()const;

    const size_t& pathToVertex(size_t vertexGlobalID) const{
        assert(vertexGlobalID>=minPathsVertex);
        assert(vertexGlobalID-minPathsVertex<vertexToPath.size());
        return vertexToPath[vertexGlobalID-minPathsVertex];
    }

    std::vector<size_t> getTimeLayersBeforePaths()const;
    std::vector<size_t> getTimeLayersAfterPaths()const;





private:
  const VertexGroups<>& vg;
  bool isFirstInterval;
  bool isLastInterval;
  std::vector<std::vector<size_t>> extractedPaths;
  size_t minPathsVertex;
  size_t maxPathsVertex;
  size_t minPathsTime;
  size_t maxPathsTime;
  std::vector<size_t> vertexToPath;
  size_t minIntervalVertex;
  //size_t lastFreeBeforePaths;  //Do not use these, they do not need to exist!
  //size_t firstFreeAfterPaths;

};


}


#endif // LDP_PATHS_EXTRACTOR_HXX
