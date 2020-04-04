#pragma once

#include <lifted_disjoint_paths/ldp_config.hxx>
#include <lifted_disjoint_paths/ldp_instance.hxx>
#include <lifted_disjoint_paths/ldp_functions.hxx>

namespace LPMP {

    namespace lifted_disjoint_paths {

    LdpInstance parse_file(const std::string& filename);
       // LdpInstance parse_file(const std::string& filename);
        //lifted_disjoint_paths_instance parse_file(const std::string& filename);
        //lifted_disjoint_paths_instance parse_string(const std::string& filename);

    } 
}
