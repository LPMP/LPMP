#pragma once

#include "discrete_tomography_instance.h"
#include <string>

namespace LPMP {

    namespace discrete_tomography_UAI_input {

        discrete_tomography_instance parse_file(const std::string& filename);
        discrete_tomography_instance parse_string(const std::string& input);

    }
}