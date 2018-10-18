#include "graph_matching/matching_problem_input.h"
#include "graph_matching/graph_matching_input.h"
#include "pegtl_parse_rules.h"
#include <cassert>

namespace LPMP {

namespace Torresani_et_al_multigraph_matching_input {

    // file consists of multiple graph matching problems. Each graph matching section beings with
    // gm x y
    // where x and y are the graph numbers (0,...)
    // then comes the graph matching problem in Torresanit et al's format.

   using parsing::opt_whitespace;
   using parsing::mand_whitespace;
   using parsing::positive_integer;

    struct graph_matching_line : pegtl::seq< opt_whitespace, pegtl::string<'g','m'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace, pegtl::eol > {};
    struct graph_matching : pegtl::star<pegtl::not_at<graph_matching_line>, pegtl::any> {};

    struct grammar :
        pegtl::star<
            pegtl::star<opt_whitespace, pegtl::eol>, 
            graph_matching_line,
            graph_matching
        >
    {};

   template< typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<> struct action< graph_matching_line > {
      template<typename INPUT>
      static void apply(const INPUT& in, multigraph_matching_input& input)
      {
         std::istringstream iss(in.string());
         char l; 
         iss >> l; assert(l == 'g');
         iss >> l; assert(l == 'm');
         std::size_t left_graph_no; iss >> left_graph_no;
         std::size_t right_graph_no; iss >> right_graph_no;
         assert(left_graph_no < right_graph_no);

         input.push_back({});
         input.back().left_graph_no = left_graph_no;
         input.back().right_graph_no = right_graph_no;
      }
   };

   template<> struct action< graph_matching > {
      template<typename INPUT>
      static void apply(const INPUT& in, multigraph_matching_input& input)
      {
          graph_matching_input gm_input = TorresaniEtAlInput::parse_string(in.string());
          input.back().gm_input = std::move(gm_input);
      }
   };

   multigraph_matching_input parse_file(const std::string& filename)
   {
      multigraph_matching_input input;
      pegtl::file_parser problem(filename);
      std::cout << "parsing " << filename << "\n";

      const bool read_success = problem.parse< grammar, action >( input );
      if(!read_success) throw std::runtime_error(std::string("could not read multigraph matching problem from file ") + filename);

      return input;
   }

   multigraph_matching_input parse_string(const std::string& problem_input)
   {
      multigraph_matching_input input;
      const bool read_success = pegtl::parse<grammar, action>(problem_input,"",input);
      if(!read_success) throw std::runtime_error("could not read multigraph matching problem from string");
      return input;
   }

} // namespace Torresani_et_al_multigraph_matching_input

} // namespace LPMP
