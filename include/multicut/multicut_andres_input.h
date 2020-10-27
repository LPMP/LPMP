#ifndef LPMP_MULTICUT_ANDRES_INPUT_H
#define LPMP_MULTICUT_ANDRES_INPUT_H

#include <string>
#include "multicut_instance.h"

namespace LPMP {

namespace multicut_andres_input {

   multicut_instance parse_file(const std::string& filename);
   multicut_instance parse_file_2d_grid_graph(const std::string& filename);

}

}

#endif // LPMP_MULTICUT_ANDRES_INPUT_H
