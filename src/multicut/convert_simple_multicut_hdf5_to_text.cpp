#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cassert>
#include <H5Cpp.h>
#include "multicut/multicut_instance.h"

using namespace std;
using namespace H5;

int main(int argc, char** argv)
{
   if(argc != 3)
      throw std::runtime_error("two input arguments expected: input file, output file.");

   const std::string input_file = argv[1];

   // Open HDF5 file handle, read only
   H5File fp(input_file.c_str(),H5F_ACC_RDONLY);

   const string edge_values_path = "/edge-values";
   const string edge_ids_path = "/graph/edges";

   // access the required dataset by path name
   DataSet edge_values_dset = fp.openDataSet(edge_values_path.c_str());
   DataSet edge_ids_dset = fp.openDataSet(edge_ids_path.c_str());

   // get the dataspace
   DataSpace edge_values_dspace = edge_values_dset.getSpace();
   DataSpace edge_ids_dspace = edge_ids_dset.getSpace();

   // get the size of the dataset
   hsize_t edge_values_dims[2];
   const hsize_t edge_values_rank = edge_values_dspace.getSimpleExtentDims(edge_values_dims, NULL); // rank = 1
   assert(edge_values_rank == 1);

   hsize_t edge_ids_dims[3];
   const hsize_t edge_ids_rank = edge_ids_dspace.getSimpleExtentDims(edge_ids_dims, NULL); // rank = 2
   assert(edge_ids_rank == 2);
   assert(edge_ids_dims[0] == edge_values_dims[0]);
   assert(edge_ids_dims[1] == 2);

   // Define the memory dataspace
   hsize_t edge_values_dimsm[1];
   edge_values_dimsm[0] = edge_values_dims[0];
   DataSpace edge_values_memspace (1,edge_values_dimsm);

   hsize_t edge_ids_dimsm[1];
   edge_ids_dimsm[0] = edge_ids_dims[0];
   edge_ids_dimsm[1] = edge_ids_dims[1];
   DataSpace edge_ids_memspace (2,edge_ids_dimsm);

   // create a vector the same size as the dataset
   vector<double> edge_values(edge_values_dims[0]);
   vector<array<std::size_t,2>> edge_ids(edge_ids_dims[0]);

   edge_values_dset.read(edge_values.data(), PredType::NATIVE_DOUBLE, edge_values_memspace, edge_values_dspace);
   edge_ids_dset.read(edge_ids.data(), PredType::STD_U64LE, edge_ids_memspace, edge_ids_dspace);

   // close the HDF5 file
   fp.close();


   LPMP::multicut_instance output;
   for(std::size_t e=0; e<edge_ids.size(); ++e)
      output.add_edge(edge_ids[e][0], edge_ids[e][1], edge_values[e]);

   const std::string output_file = argv[2];
   std::ofstream file_stream(output_file, std::ofstream::out);
   output.write_problem(file_stream);
   file_stream.close();
}
