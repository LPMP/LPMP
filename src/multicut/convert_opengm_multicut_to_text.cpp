#include <fstream>
#include <vector>
#include <string>
#include "opengm/opengm.hxx"
#include "opengm/graphicalmodel/graphicalmodel.hxx"
#include "opengm/graphicalmodel/graphicalmodel_hdf5.hxx"
#include "opengm/operations/adder.hxx"
#include "assert.h"

std::vector<std::string> split(const std::string &s, const char delim) {
   std::vector<std::string> elems;
   std::stringstream ss(s);
   std::string item;
   while (std::getline(ss, item, delim)) {
      if (item.length() > 0) {
         elems.push_back(item);  
      }
   }
   return elems;
}


int main(int argc, char** argv)
{
   assert(argc == 3); // first arg input in opengm format, second is output in text format

   auto input = split(argv[1],':');
   if(input.size() > 2) throw std::runtime_error("invalid input.");
   const std::string input_file = input[1];
   std::string input_model;
   if(input.size() == 1) {
      input_model = "gm";
   } else {
      input_model = input[1];
   }

   std::cout << "open file " << input_file << " with model " << input_model << "\n";
   
   const std::string output_file = argv[2];

   typedef double                                                               ValueType;          // type used for values
   typedef size_t                                                               IndexType;          // type used for indexing nodes and factors (default : size_t)
   typedef size_t                                                               LabelType;          // type used for labels (default : size_t)
   typedef opengm::Adder                                                        OpType;             // operation used to combine terms
   typedef opengm::ExplicitFunction<ValueType,IndexType,LabelType>              ExplicitFunction;   // shortcut for explicite function
   typedef opengm::PottsFunction<ValueType,IndexType,LabelType>                 PottsFunction;      // Potts function
   typedef opengm::meta::TypeListGenerator<ExplicitFunction,PottsFunction>::type              FunctionTypeList;   // list of all function the model cal use (this trick avoids virtual methods) - here only one
   typedef opengm::DiscreteSpace<IndexType, LabelType>                          SpaceType;          // type used to define the feasible statespace
   typedef opengm::GraphicalModel<ValueType,OpType,FunctionTypeList,SpaceType>  Model;              // type of the model
   typedef Model::FunctionIdentifier                                            FunctionIdentifier; // type of the function identifier

   Model gm;
   opengm::hdf5::load(gm, input_file,input_model);

   std::ofstream file_stream(output_file, std::ofstream::out);
   file_stream << "MULTICUT\n";

   for(std::size_t e=0; e<gm.numberOfFactors(); ++e) {
     assert(gm[e].numberOfVariables() == 2);
     const std::size_t i = gm.variableOfFactor(e,0);
     const std::size_t j = gm.variableOfFactor(e,1);
     assert(i < j);

     const double cost = gm[e](std::array<std::size_t,2>({1,0}).begin()) - gm[e](std::array<std::size_t,2>({0,0}).begin());

     file_stream << i << " " << j << " " << cost << "\n";
   }

   file_stream.close();
}
