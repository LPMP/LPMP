#include "max_cut/max_cut_text_input.h"
#include "cut_base/cut_base_text_input.hxx"

namespace LPMP {

    namespace max_cut_text_input {

        using max_cut_identifier = pegtl::opt<pegtl::string<'M','A','X','-','C','U','T'>>;

        max_cut_instance parse_string(const std::string& input)
        {
            return cut_base_text_input::parse_string<max_cut_instance, max_cut_identifier, -1>(input);
        }
        max_cut_instance parse_file(const std::string& filename)
        {
            return cut_base_text_input::parse_file<max_cut_instance, max_cut_identifier, -1>(filename);
        }

    }

}
