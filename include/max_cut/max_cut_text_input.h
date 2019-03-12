#pragma once

#include "max_cut_instance.hxx"

namespace LPMP {
    namespace max_cut_text_input {
        max_cut_instance parse_string(const std::string& input);
        max_cut_instance parse_file(const std::string& filename); 
    }
}
