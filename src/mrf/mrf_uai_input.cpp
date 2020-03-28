#include "mrf/mrf_uai_input.h"
#include "mrf/mrf_uai_input_impl.hxx"

namespace LPMP {

namespace binary_MRF_uai_input {

   binary_MRF_instance parse_file(const std::string& filename)
   {
       return mrf_uai_input::parse_file<pegtl::string<'M','A','R','K','O','V'>, binary_MRF_instance>(filename); 
   }
   binary_MRF_instance parse_string(const std::string& uai_string)
   {
       return mrf_uai_input::parse_string<pegtl::string<'M','A','R','K','O','V'>, binary_MRF_instance>(uai_string);
   }

}

namespace mrf_uai_input {

   mrf_input parse_file(const std::string& filename)
   {
       return parse_file<pegtl::string<'M','A','R','K','O','V'>>(filename);
   }

   mrf_input parse_string(const std::string& uai_string)
   {
       return parse_string<pegtl::string<'M','A','R','K','O','V'>>(uai_string);
   }


} // namespace mrf_uai_input

} // namespace LPMP
