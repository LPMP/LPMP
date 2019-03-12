#ifndef LPMP_MULTICUT_OPENGM_INPUT
#define LPMP_MULTICUT_OPENGM_INPUT

#include <string>
#include "multicut_instance.h"

namespace LPMP {

namespace multicut_opengm_input {
   multicut_instance parse_file(const std::string& filename);
}

}

#endif // LPMP_MULTICUT_OPENGM_INPUT

