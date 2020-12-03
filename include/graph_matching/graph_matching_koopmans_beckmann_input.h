#pragma once

#include <vector>
#include <string>
#include <iostream>
#include "graph_matching_instance_koopmans_beckmann.h"

namespace LPMP {

namespace graph_matching_koopmans_beckmann_input {

   graph_matching_instance_koopmans_beckmann parse_file(const std::string& filename);

   graph_matching_instance_koopmans_beckmann parse_string(const std::string& input);

}

} // namespace LPMP

