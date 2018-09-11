#include "graph_matching/multigraph_matching_input.h"
#include "pegtl_parse_rules.h"
#include <cassert>

namespace LPMP {

namespace multigraph_matching_input {

    // file consists of multiple graph matching problems. Each graph matching section beings with
    // gm x y
    // where x and y are the graph numbers (0,...)
    // then comes the graph matching problem in Torresanit et al's format.

   using parsing::opt_whitespace;
   using parsing::mand_whitespace;
   using parsing::positive_integer;

    struct graph_matching_line : pegtl::seq< opt_whitespace, pegtl::string<'g','m'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace, pegtl::eol > {};
    struct graph_matching : pegtl::star<pegtl::not_at<graph_matching_line>, pegtl::any> {};
    //: pegtl::until<graph_matching_line> {};

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
      static void apply(const INPUT& in, mgm_input& input)
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
      static void apply(const INPUT& in, mgm_input& input)
      {
          graph_matching_input gm_input = TorresaniEtAlInput::parse_string(in.string());
          input.back().gm_input = gm_input;
      }
   };

   mgm_input parse_file(const std::string& filename)
   {
      mgm_input input;
      pegtl::file_parser problem(filename);
      std::cout << "parsing " << filename << "\n";

      const bool read_success = problem.parse< grammar, action >( input );
      assert(read_success);

      return input;
   }

   mgm_input parse_string(const std::string& problem_input)
   {
      mgm_input input;
      const bool read_success = pegtl::parse<grammar, action>(problem_input,"",input);
      assert(read_success);
      return input;
   }

} // namespace multigraph_matching_input 

} // namespace LPMP
