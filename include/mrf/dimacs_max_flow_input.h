#ifndef LPMP_DIMACS_MAX_FLOW_INPUT_H
#define LPMP_DIMACS_MAX_FLOW_INPUT_H

#include <string>
#include "max_flow_instance.hxx"

namespace LPMP {

    namespace dimacs_max_flow_input {

        max_flow_instance parse_file(const std::string& filename);
        max_flow_instance parse_string(const std::string& string);

    }

}

#endif // LPMP_DIMACS_MAX_FLOW_INPUT_H
