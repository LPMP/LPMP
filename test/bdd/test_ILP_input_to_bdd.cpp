#include "test.h"
#include "cuddObj.hh"
#include "bdd/convert_pb_to_bdd.h"
#include <vector>

using namespace LPMP;

void write_dd(DdManager *gbm, DdNode *dd, const char* filename)
{
    FILE *outfile; // output file pointer for .dot file
    outfile = fopen(filename,"w");
    DdNode **ddnodearray = (DdNode**)malloc(sizeof(DdNode*)); // initialize the function array
    ddnodearray[0] = dd;
    Cudd_DumpDot(gbm, 1, ddnodearray, NULL, NULL, outfile); // dump the function to .dot file
    free(ddnodearray);
    fclose (outfile); // close the file */
}

std::vector<int> create_vector(const std::vector<std::size_t> one_indices, const std::vector<std::size_t> minus_one_indices)
{
    std::vector<int> vec;
    for(auto x : one_indices) {
        if(x >= vec.size())
            vec.resize(x+1, 0.0);
        assert(vec[x] == 0.0);
        vec[x] = 1.0;
    }
    for(auto x : minus_one_indices) {
        if(x >= vec.size())
            vec.resize(x+1, 0.0);
        assert(vec[x] == 0.0);
        vec[x] = -1.0;
    }

    return vec;
}

int main(int argc, char** argv)
{
    Cudd bdd_mgr;
    bdd_converter converter(bdd_mgr);

    {
        std::vector<int> simplex_weights = {1,1,1,1,1};
        auto simplex_bdd = converter.convert_to_bdd(simplex_weights.begin(), simplex_weights.end(), inequality_type::equal, 1);

        std::vector<BDD> bdds = {simplex_bdd};

        //auto add = Cudd_BddToAdd(bdd_mgr.getManager(), simplex_bdd.getNode());
        //write_dd(bdd_mgr.getManager(), add, "kwas.dot");
        bdd_mgr.DumpDot( bdds );
        //char const * const * inames = 0, 
        //      char const * const * onames = 0, 
        //      FILE * fp = stdout) const;
    }

    // mrf with three variables and Potts potentials
    { 
        return 0;
        Cudd bdd_mgr;
        auto simplex_1 = converter.convert_to_bdd(create_vector({0,1,2,3,4},{}), inequality_type::equal, 1);
        auto simplex_2 = converter.convert_to_bdd(create_vector({5,6,7,8,9},{}), inequality_type::equal, 1);
        auto simplex_3 = converter.convert_to_bdd(create_vector({10,11,12,13,14},{}), inequality_type::equal, 1);

        auto potts_12_1 = converter.convert_to_bdd(create_vector({0},{5,15}), inequality_type::smaller_equal , 0);
        auto potts_12_2 = converter.convert_to_bdd(create_vector({1},{6,15}), inequality_type::smaller_equal , 0);
        auto potts_12_3 = converter.convert_to_bdd(create_vector({2},{7,15}), inequality_type::smaller_equal , 0);
        auto potts_12_4 = converter.convert_to_bdd(create_vector({3},{8,15}), inequality_type::smaller_equal , 0);
        auto potts_12_5 = converter.convert_to_bdd(create_vector({4},{9,15}), inequality_type::smaller_equal , 0);
        auto potts_21_1 = converter.convert_to_bdd(create_vector({5},{0,15}), inequality_type::smaller_equal , 0);
        auto potts_21_2 = converter.convert_to_bdd(create_vector({6},{1,15}), inequality_type::smaller_equal , 0);
        auto potts_21_3 = converter.convert_to_bdd(create_vector({7},{2,15}), inequality_type::smaller_equal , 0);
        auto potts_21_4 = converter.convert_to_bdd(create_vector({8},{3,15}), inequality_type::smaller_equal , 0);
        auto potts_21_5 = converter.convert_to_bdd(create_vector({9},{4,15}), inequality_type::smaller_equal , 0);

        auto potts_13_1 = converter.convert_to_bdd(create_vector({0},{10,16}), inequality_type::smaller_equal , 0);
        auto potts_13_2 = converter.convert_to_bdd(create_vector({1},{11,16}), inequality_type::smaller_equal , 0);
        auto potts_13_3 = converter.convert_to_bdd(create_vector({2},{12,16}), inequality_type::smaller_equal , 0);
        auto potts_13_4 = converter.convert_to_bdd(create_vector({3},{13,16}), inequality_type::smaller_equal , 0);
        auto potts_13_5 = converter.convert_to_bdd(create_vector({4},{14,16}), inequality_type::smaller_equal , 0);
        auto potts_31_1 = converter.convert_to_bdd(create_vector({10},{0,16}), inequality_type::smaller_equal , 0);
        auto potts_31_2 = converter.convert_to_bdd(create_vector({11},{1,16}), inequality_type::smaller_equal , 0);
        auto potts_31_3 = converter.convert_to_bdd(create_vector({12},{2,16}), inequality_type::smaller_equal , 0);
        auto potts_31_4 = converter.convert_to_bdd(create_vector({13},{3,16}), inequality_type::smaller_equal , 0);
        auto potts_31_5 = converter.convert_to_bdd(create_vector({14},{4,16}), inequality_type::smaller_equal , 0);

        auto potts_23_1 = converter.convert_to_bdd(create_vector({5},{10,17}), inequality_type::smaller_equal , 0);
        auto potts_23_2 = converter.convert_to_bdd(create_vector({6},{11,17}), inequality_type::smaller_equal , 0);
        auto potts_23_3 = converter.convert_to_bdd(create_vector({7},{12,17}), inequality_type::smaller_equal , 0);
        auto potts_23_4 = converter.convert_to_bdd(create_vector({8},{13,17}), inequality_type::smaller_equal , 0);
        auto potts_23_5 = converter.convert_to_bdd(create_vector({9},{14,17}), inequality_type::smaller_equal , 0);
        auto potts_32_1 = converter.convert_to_bdd(create_vector({10},{5,17}), inequality_type::smaller_equal , 0);
        auto potts_32_2 = converter.convert_to_bdd(create_vector({11},{6,17}), inequality_type::smaller_equal , 0);
        auto potts_32_3 = converter.convert_to_bdd(create_vector({12},{7,17}), inequality_type::smaller_equal , 0);
        auto potts_32_4 = converter.convert_to_bdd(create_vector({13},{8,17}), inequality_type::smaller_equal , 0);
        auto potts_32_5 = converter.convert_to_bdd(create_vector({14},{9,17}), inequality_type::smaller_equal , 0);

        auto all_simplex = simplex_1.And(simplex_2.And( simplex_3));
        auto bdd_potts_12 = potts_12_1.And(potts_12_2.And(potts_12_3.And(potts_12_4.And(potts_12_5.And(potts_21_1.And(potts_21_2.And(potts_21_3.And(potts_21_4.And(potts_21_5)))))))));
        auto bdd_potts_13 = potts_13_1.And(potts_13_2.And(potts_13_3.And(potts_13_4.And(potts_13_5.And(potts_31_1.And(potts_31_2.And(potts_31_3.And(potts_31_4.And(potts_31_5)))))))));
        auto bdd_potts_23 = potts_23_1.And(potts_23_2.And(potts_23_3.And(potts_23_4.And(potts_23_5.And(potts_32_1.And(potts_32_2.And(potts_32_3.And(potts_32_4.And(potts_32_5)))))))));

        auto all = all_simplex.And(bdd_potts_12.And(bdd_potts_13.And(bdd_potts_23)));

        bdd_mgr.ReduceHeap(CUDD_REORDER_EXACT);

        std::vector<BDD> bdds = {all};
        bdd_mgr.DumpDot( bdds ); 
    }
}
