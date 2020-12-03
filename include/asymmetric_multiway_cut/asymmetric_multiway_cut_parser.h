#pragma once

#include "asymmetric_multiway_cut_instance.h"

namespace LPMP {

    namespace asymmetric_multiway_cut_parser {

        asymmetric_multiway_cut_instance parse_file(const std::string& filename);
        asymmetric_multiway_cut_instance parse_string(const std::string& string);

    }
}
