#ifndef LPMP_MRF_UAI_INPUT_IMPL_HXX
#define LPMP_MRF_UAI_INPUT_IMPL_HXX

#include "mrf_input.h"
#include "pegtl_parse_rules.h"
#include "pegtl/parse.hh"
#include <variant>

namespace LPMP {

namespace mrf_uai_input {

   // import basic parsers
   using parsing::opt_whitespace;
   using parsing::mand_whitespace;
   using parsing::opt_invisible;
   using parsing::mand_invisible;
   using parsing::positive_integer;
   using parsing::real_number;

   struct mrf_input_helper {
       std::vector<std::size_t> cardinality;
       using clique_scope_type = std::variant< std::size_t, std::array<std::size_t,2>, std::vector<std::size_t> >;
       std::vector<clique_scope_type> clique_scopes;
       std::size_t number_of_cliques;

       std::vector<double> current_function_table;
       std::size_t current_function_table_size;
       std::size_t current_clique_number = 0;
       std::size_t current_pairwise_clique_number = 0;
   };

   struct number_of_variables : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace > {};
   // vector of integers denoting how many labels each variable has
   struct cardinality : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace > {};
   struct allocate_data : pegtl::seq<> {};
   struct number_of_cliques : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace> {};
   // first is the number of variables in the clique, then the actual variables.
   // the clique_scopes should match number_of_clique_lines, each line consisting of a sequence of integers
   struct new_clique_scope : pegtl::seq< positive_integer > {};
   struct clique_scope : pegtl::seq< positive_integer > {};
   struct clique_scope_line : pegtl::seq< opt_whitespace, new_clique_scope, pegtl::plus< opt_whitespace, clique_scope >, opt_whitespace, pegtl::eol > {};
   struct clique_scopes_end
   {
      template< pegtl::apply_mode A, template< typename ... > class Action, template< typename ... > class Control, typename Input >
         static bool match(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
         {
            return input_helper.number_of_cliques == input_helper.clique_scopes.size();
         }
   };
   struct clique_scopes : pegtl::until< clique_scopes_end, clique_scope_line > {};
   // a function table is begun by number of entries and then a list of real numbers. Here we record all the values in the real stack
   // do zrobienia: treat whitespace

   struct new_function_table : pegtl::seq< positive_integer > {};
   struct function_table_entry : pegtl::seq< real_number > {};
   struct function_tables_end
   {
      template< pegtl::apply_mode A, template< typename ... > class Action, template< typename ... > class Control, typename Input >
         static bool match(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
         {
             return input_helper.number_of_cliques == input_helper.current_clique_number;
         }

   };
   struct function_table_end
   {
      template< pegtl::apply_mode A, template< typename ... > class Action, template< typename ... > class Control, typename Input >
         static bool match(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
         {
             if(input_helper.current_function_table.size() == input_helper.current_function_table_size) {

                auto& table = input_helper.current_function_table;
                const std::size_t clique_number = input_helper.current_clique_number;
                auto& clique_scope = input_helper.clique_scopes[clique_number];
                // write function table into unary, pairwise or higher order table
                if( std::holds_alternative<std::size_t>( clique_scope ) ) {
                    const auto var = std::get<std::size_t>(clique_scope);
                    assert(input.cardinality(var) == table.size());
                    for(std::size_t i=0; i<table.size(); ++i) {
                        input.unaries[var][i] = table[i];
                    }
                } else if( std::holds_alternative<std::array<std::size_t,2>>( clique_scope ) ) {
                    auto& vars = std::get<std::array<std::size_t,2>>(clique_scope);
                    assert(input.cardinality(vars[0])*input.cardinality(vars[1]) == table.size());
                    for(std::size_t i=0; i<input.cardinality(vars[0]); ++i) {
                        for(std::size_t j=0; j<input.cardinality(vars[1]); ++j) {
                            input.pairwise_values(input_helper.current_pairwise_clique_number, i, j) = table[i*input.cardinality(vars[1]) + j];
                        }
                    }
                    input_helper.current_pairwise_clique_number++; 
                } else {
                    assert(std::holds_alternative<std::vector<std::size_t>>(clique_scope));
                    assert(false); // not implemented yet
                } 

                table.clear();
                input_helper.current_clique_number++;
                return true;

            } else {
                return false;
            }
         }
   };
   struct function_table : pegtl::seq< new_function_table, opt_invisible, pegtl::until< function_table_end, opt_invisible, function_table_entry >, opt_invisible > {};
   struct function_tables : pegtl::seq< opt_invisible, pegtl::until< function_tables_end, function_table >, opt_invisible > {};

   template< typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<> struct action< number_of_variables > {
      template<typename Input>
      static void apply(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
      {
         input_helper.cardinality.reserve(std::stoul(in.string()));
      }
   };

   template<> struct action< number_of_cliques > {
      template<typename Input>
      static void apply(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
      {
         input_helper.number_of_cliques = std::stoul(in.string()); 
      }
   };

   template<> struct action< cardinality > {
      template<typename Input>
      static void apply(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
      {
         input_helper.cardinality.push_back(std::stoul(in.string()));
      }
   };

   template<> struct action< allocate_data > {
      template<typename Input>
      static void apply(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
      {
          // allocate space for unary and pairwise factors
          input.unaries.resize(input_helper.cardinality.begin(), input_helper.cardinality.end());
          // we must write 0 into unaries, since not all unary potentials must be present
          for(std::size_t i=0; i<input.unaries.size(); ++i) {
              for(std::size_t j=0; j<input.unaries[i].size(); ++j) {
                  input.unaries(i,j) = 0.0;
              }
          }

          std::vector<std::array<std::size_t,2>> pairwise_size;
          pairwise_size.reserve(input.pairwise_indices.size());
          for(const auto pairwise_index : input.pairwise_indices) {
              pairwise_size.push_back( {input.cardinality(pairwise_index[0]), input.cardinality(pairwise_index[1])} );
          }

          input.pairwise_values.resize(pairwise_size.begin(), pairwise_size.end());
      } 
   };

   template<> struct action< new_clique_scope > {
      template<typename Input>
      static void apply(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
      {
          const std::size_t arity = std::stoul(in.string());
          assert(arity == 1 || arity == 2); // higher order not yet supported
          mrf_input_helper::clique_scope_type clique_scope;
          if(arity == 1) {
              clique_scope = std::size_t(std::numeric_limits<std::size_t>::max());
          } else if(arity == 2) {
              clique_scope = std::array<std::size_t,2>( {std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max()} );
          } else {
              assert(false);
          }
          input_helper.clique_scopes.push_back( clique_scope );
      }
   };

   template<> struct action< clique_scope > {
      template<typename Input>
      static void apply(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
      {
          const auto var_no = std::stoul(in.string());
          assert(var_no < input_helper.cardinality.size());
          auto& clique_scope = input_helper.clique_scopes.back();
          if( std::holds_alternative<std::size_t>( clique_scope ) ) {
              std::get<std::size_t>(clique_scope) = var_no;
          } else if( std::holds_alternative<std::array<std::size_t,2>>( clique_scope ) ) {
              auto& pairwise_idx = std::get<std::array<std::size_t,2>>(clique_scope);
              if(pairwise_idx[0] == std::numeric_limits<std::size_t>::max()) {
                  pairwise_idx[0] = var_no;
              } else {
                  assert(pairwise_idx[1] == std::numeric_limits<std::size_t>::max());
                  pairwise_idx[1] = var_no;
                  input.pairwise_indices.push_back({pairwise_idx[0], pairwise_idx[1]});
              }
          } else {
              assert(std::holds_alternative<std::vector<std::size_t>>(clique_scope));
              std::get<std::vector<std::size_t>>(clique_scope).push_back(var_no);
          }
      }
   };
   template<> struct action< new_function_table > {
      template<typename Input>
      static void apply(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
      {
         input_helper.current_function_table.clear();
         input_helper.current_function_table_size  = std::stoul(in.string());
         input_helper.current_function_table.reserve( input_helper.current_function_table_size );
      }
   };
   template<> struct action< function_table_entry > {
      template<typename Input>
      static void apply(const Input& in, mrf_input& input, mrf_input_helper& input_helper) 
      {
         input_helper.current_function_table.push_back(std::stod(in.string()));
      }
   };

   template<typename PEGTL_STRING>
   struct grammar {
       struct type : pegtl::seq<
                        opt_whitespace, PEGTL_STRING, opt_whitespace, pegtl::eol,
                        number_of_variables, pegtl::eol,
                        pegtl::plus< cardinality >, pegtl::eol,
                        number_of_cliques, pegtl::eol,
                        clique_scopes,
                        opt_invisible,
                        allocate_data,
                        function_tables
                        > {};
   }; 

   template<typename PEGTL_STRING>
   mrf_input parse_file(const std::string& filename)
   {
       mrf_input input;
       mrf_input_helper input_helper;
       pegtl::file_parser problem(filename);

       const bool read_success = problem.parse<typename grammar<PEGTL_STRING>::type, action>( input, input_helper );
       assert(read_success);
       return input; 
   }

   template<typename PEGTL_STRING>
   mrf_input parse_string(const std::string& uai_string)
   {
       mrf_input input;
       mrf_input_helper input_helper;

       const bool read_success = pegtl::parse<typename grammar<PEGTL_STRING>::type, action>(uai_string, "", input, input_helper);
       assert(read_success);
       return input;
   }

} // namespace mrf_uai_input

} // namespace LPMP

#endif // LPMP_MRF_UAI_INPUT_IMPL_HXX

