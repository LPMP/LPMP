#ifndef LPMP_MRF_OPENGM_INPUT_H
#define LPMP_MRF_OPENGM_INPUT_H

#include "mrf_input.h"

#include <string>
#include <cassert>

namespace LPMP {

// file format described in http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php
namespace mrf_opengm_input {

   mrf_input parse_file(const std::string& filename);
   mrf_input parse_string(const std::string& filename);

} // namespace mrf_opengm_input

} // namespace LPMP

#endif // LPMP_MRF_OPENGM_INPUT_H

