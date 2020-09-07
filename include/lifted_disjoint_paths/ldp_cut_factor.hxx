#ifndef LDP_CUT_FACTOR_HXX
#define LDP_CUT_FACTOR_HXX

#endif // LDP_CUT_FACTOR_HXX
#include <utility>
#include<cstdlib>
#include<vector>
#include<array>
#include<set>

namespace LPMP {


class ldp_cut_factor
{
public:
    ldp_cut_factor(size_t v_,size_t w_, std::vector<std::array<size_t,2>> inputEdges,std::vector<double> inputCosts) { //TODO: maybe inputEdges as map<size_t<map<size_t,double>> - better for creating edges as pairs of order of v1 and v2
        std::set<size_t> inpVertices;
        std::set<size_t> outVertices;
        for (size_t i = 0; i < inputEdges.size(); ++i) {
            size_t v1=inputEdges[i][0];
            size_t v2=inputEdges[i][1];
            if(inpVertices.count(v1)==0){
                inpVertices.insert(v1);
            }

            if(outVertices.count(v2)==0){
                outVertices.insert(v2);
            }

        }
        //inputVertices=std::vector<size_t>(inpVertices.size());
        for(size_t vert:inpVertices){
            inputVertices.push_back(vert);
        }
        for(size_t vert:outVertices){
            outputVertices.push_back(vert);
        }

    }


private:
std::vector<double> costs;
std::vector<std::array<size_t,2>> edges;
std::vector<size_t> inputVertices;
std::vector<size_t> outputVertices;
size_t v;
size_t w;





};



}
