#pragma once

#include <string>
#include "multicut_instance.h"

namespace LPMP {

   namespace multicut_text_input {

      multicut_instance parse_string(const std::string& input);
      multicut_instance parse_file(const std::string& input);

   }

}
