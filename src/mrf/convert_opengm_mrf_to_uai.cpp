#include <fstream>
#include "hdf5.h"
#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>
#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/functions/explicit_function.hxx>
#include <opengm/functions/potts.hxx>
#include <opengm/functions/pottsn.hxx>
#include <opengm/functions/pottsg.hxx>
#include "opengm/functions/truncated_absolute_difference.hxx"
#include "opengm/functions/truncated_squared_difference.hxx"

int main(int argc, char** argv)
{
   assert(argc == 3); // first arg input in opengm format, second is output in text format
   const std::string input_file = argv[1];
   const std::string output_file = argv[2];

   typedef double ValueType;
   typedef size_t IndexType;
   typedef size_t LabelType;
   typedef opengm::Adder OperatorType;
   typedef opengm::Minimizer AccumulatorType;
   typedef opengm::DiscreteSpace<IndexType, LabelType> SpaceType;

   // Set functions for graphical model
   typedef opengm::meta::TypeListGenerator<
      opengm::ExplicitFunction<ValueType, IndexType, LabelType>,
      opengm::PottsFunction<ValueType, IndexType, LabelType>,
      opengm::PottsNFunction<ValueType, IndexType, LabelType>,
      opengm::PottsGFunction<ValueType, IndexType, LabelType>,
      opengm::TruncatedSquaredDifferenceFunction<ValueType, IndexType, LabelType>,
      opengm::TruncatedAbsoluteDifferenceFunction<ValueType, IndexType, LabelType>
   >::type FunctionTypeList;


   typedef opengm::GraphicalModel<
      ValueType,
      OperatorType,
      FunctionTypeList,
      SpaceType
   > GmType;
   

   GmType gm; 
   opengm::hdf5::load(gm, input_file,"gm");

   std::ofstream uai(output_file, std::ofstream::out);
   uai << "MARKOV\n";
   uai << gm.numberOfVariables() << "\n";
   for(std::size_t i=0; i<gm.numberOfVariables(); ++i) {
     uai << gm.numberOfLabels(i) << " ";
   }
   uai << "\n";
   uai << gm.numberOfFactors() << "\n";

   for(std::size_t f=0; f<gm.numberOfFactors(); ++f){
     if(!(gm[f].numberOfVariables() == 1 || gm[f].numberOfVariables() == 2)) {
       std::cout << "graphical models with unary and pairwise variables only supported\n";
       exit(-1);
     }
     if(gm[f].numberOfVariables()==1){
       uai << "1 " << gm.variableOfFactor(f,0) << "\n";
     }
     if(gm[f].numberOfVariables()==2){
       const std::size_t i = gm.variableOfFactor(f,0);
       const std::size_t j = gm.variableOfFactor(f,1);
       uai << "2 " << i << " " << j << "\n";
     }
   }

   for(std::size_t f=0; f<gm.numberOfFactors(); ++f){
     if(gm[f].numberOfVariables()==1){
       const std::size_t i = gm.variableOfFactor(f,0);
       uai << gm.numberOfLabels(i) << "\n";
       for(std::size_t l=0; l<gm[f].numberOfLabels(0); ++l){
         uai << gm[f](std::array<std::size_t,1>({l}).begin()) << " ";
       } 
       uai << "\n\n";
     }
     if(gm[f].numberOfVariables()==2){
       const std::size_t i = gm.variableOfFactor(f,0);
       const std::size_t j = gm.variableOfFactor(f,1);
       uai << gm.numberOfLabels(i)*gm.numberOfLabels(j) << "\n";
       for(std::size_t l1=0; l1<gm[f].numberOfLabels(0); ++l1){
         for(std::size_t l2=0; l2<gm[f].numberOfLabels(1); ++l2){
           uai << gm[f](std::array<std::size_t,2>({l1,l2}).begin()) << " ";
         }
         uai << "\n";
       }
       uai << "\n";
     }
   }
}
