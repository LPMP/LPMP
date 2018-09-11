#ifndef LPMP_MULTIGRAPH_MATCHING_INPUT_H
#define LPMP_MULTIGRAPH_MATCHING_INPUT_H

#include "graph_matching_input.h"
#include <vector>
#include <string>

namespace LPMP {

struct mgm_input_entry {
    std::size_t left_graph_no, right_graph_no;
    graph_matching_input gm_input;
};

using mgm_input = std::vector<mgm_input_entry>;

namespace multigraph_matching_input {

    // file consists of multiple graph matching problems. Each graph matching section beings with
    // gm x y
    // where x and y are the graph numbers (0,...)
    // then comes the graph matching problem in Torresanit et al's format.


    mgm_input parse_file(const std::string& filename);
    mgm_input parse_string(const std::string& input);

} // namespace multigraph_matching_input 

} // namespace LPMP

#endif // LPMP_MULTIGRAPH_MATCHING_INPUT_H
