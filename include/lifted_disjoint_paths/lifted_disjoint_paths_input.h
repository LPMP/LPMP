#pragma once

#include <string>
#include "lifted_disjoint_paths_instance.h"

namespace LPMP {

    namespace lifted_disjoint_paths {

        lifted_disjoint_paths_instance parse_file(const std::string& filename);
        lifted_disjoint_paths_instance parse_string(const std::string& filename);

    } 
}
