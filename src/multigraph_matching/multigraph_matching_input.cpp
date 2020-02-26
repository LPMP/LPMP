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

    struct comment_line : pegtl::seq< opt_whitespace, pegtl::sor<pegtl::string<'c'>, pegtl::string<'#'>>, pegtl::until< pegtl::eol >> {};
    struct empty_line : pegtl::seq< opt_whitespace, pegtl::eol > {};
    struct ignore_line : pegtl::sor<comment_line, empty_line > {};
    struct graph_matching_line : pegtl::seq< opt_whitespace, pegtl::string<'g','m'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace, pegtl::eol > {};
    struct graph_matching : pegtl::star<pegtl::not_at<graph_matching_line>, pegtl::any> {};

    struct grammar :
        pegtl::star<
            pegtl::star<ignore_line>, 
            graph_matching_line,
            graph_matching,
            pegtl::eof
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

         std::cout << "read graph matching problem " << left_graph_no << " -> " << right_graph_no << "\n";
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
      if(!read_success) 
          throw std::runtime_error(std::string("could not read multigraph matching problem from file ") + filename);

      return input;
   }

   multigraph_matching_input parse_string(const std::string& problem_input)
   {
      multigraph_matching_input input;
      const bool read_success = pegtl::parse<grammar, action>(problem_input,"",input);
      if(!read_success) 
          throw std::runtime_error("could not read multigraph matching problem from string");
      return input;
   }

} // namespace Torresani_et_al_multigraph_matching_input

using parsing::opt_whitespace;
using parsing::mand_whitespace;
using parsing::positive_integer;

struct graph_matching_line : pegtl::seq< opt_whitespace, pegtl::string<'g','m'>, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace, pegtl::eol > {};
struct graph_matching : pegtl::star<pegtl::not_at<graph_matching_line>, pegtl::any> {};
struct comment_line : pegtl::seq< opt_whitespace, pegtl::sor<pegtl::string<'c'>, pegtl::string<'#'>>, pegtl::until< pegtl::eol >> {};
struct empty_line : pegtl::seq< opt_whitespace, pegtl::eol > {};
struct ignore_line : pegtl::sor<comment_line, empty_line > {};

struct grammar :
   pegtl::star<
   pegtl::star<ignore_line>, 
   graph_matching_line,
   graph_matching
   >
{};

template< typename Rule >
struct action
: pegtl::nothing< Rule > {};

template<> struct action< graph_matching_line > {
   template<typename INPUT>
      static void apply(const INPUT& in, multigraph_matching_input::labeling& input)
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
      static void apply(const INPUT& in, multigraph_matching_input::labeling& input)
      {
         graph_matching_input::labeling gm_input = parse_graph_matching_result_string(in.string());
         input.back().labeling = std::move(gm_input);
      }
};

multigraph_matching_input::labeling parse_multigraph_matching_result_file(const std::string& filename)
{
   multigraph_matching_input::labeling output;
   throw std::runtime_error("not implemented yet");
   return output;
}
multigraph_matching_input::labeling parse_multigraph_matching_result_string(const std::string& input)
{
   multigraph_matching_input::labeling output;
   const bool read_success = pegtl::parse<grammar, action>(input,"",output);
   if(!read_success) throw std::runtime_error("could not read multigraph matching result from string");
   return output;
}

} // namespace LPMP
