#include "mrf/mrf_opengm_input.h"
#include <stdexcept> 
#include "hdf5.h"
#include "opengm/opengm.hxx"
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

namespace LPMP {

namespace mrf_opengm_input {

   mrf_input parse_file(const std::string& filename)
   {
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
       opengm::hdf5::load(gm, filename,"gm");

       mrf_input input;

       std::vector<std::size_t> cardinality(gm.numberOfVariables());
       for(std::size_t i=0; i<gm.numberOfVariables(); ++i) {
           cardinality[i] = gm.numberOfLabels(i);
       }
       input.unaries.resize(cardinality.begin(), cardinality.end());

       std::vector<std::array<std::size_t,2>> pairwise_size;
       for(std::size_t f=0; f<gm.numberOfFactors(); ++f) {
           if(gm[f].numberOfVariables()==2) {
               pairwise_size.push_back({gm[f].numberOfLabels(0), gm[f].numberOfLabels(1)});
               input.pairwise_indices.push_back({gm.variableOfFactor(f,0), gm.variableOfFactor(f,0)});
           }
       }
       input.pairwise_values.resize(pairwise_size.begin(), pairwise_size.end());

       std::size_t pairwise_counter = 0;
       for(std::size_t f=0; f<gm.numberOfFactors(); ++f){

           if(gm[f].numberOfVariables()==0){
               // ignore for now
           }
           else if(gm[f].numberOfVariables()==1){
               const std::size_t i = gm.variableOfFactor(f,0);
               for(std::size_t l=0; l<gm[f].numberOfLabels(0); ++l){
                   input.unaries(i,l) = gm[f](std::array<std::size_t,1>({l}).begin()); 
               } 
           } 
           else if(gm[f].numberOfVariables()==2){
               const std::size_t i = gm.variableOfFactor(f,0);
               const std::size_t j = gm.variableOfFactor(f,1);
               for(std::size_t l1=0; l1<gm[f].numberOfLabels(0); ++l1){
                   for(std::size_t l2=0; l2<gm[f].numberOfLabels(1); ++l2){
                       input.pairwise_values(pairwise_counter, l1, l2) = gm[f](std::array<std::size_t,2>({l1,l2}).begin()); 
                   }
               }
               pairwise_counter++;
           }
           else{
               throw std::runtime_error("Factors of order higher than 2 are so far not supported.");
           }

       }

       return input;
   }

} // namespace mrf_opengm_input

} // namespace LPMP
