#ifndef LPMP_MULTICUT_TEXT_INPUT
#define LPMP_MULTICUT_TEXT_INPUT

#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <cassert>
#include "multicut_instance.hxx"

namespace LPMP {

   namespace multicut_text_input {

      multicut_instance parse_string(const std::string& input);
      multicut_instance parse_file(const std::string& input);

   }

}

#endif // LPMP_MULTICUT_TEXT_INPUT
