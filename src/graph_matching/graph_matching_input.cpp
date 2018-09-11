#include "graph_matching/graph_matching_input.h"
#include "pegtl_parse_rules.h"
#include <cassert>

namespace LPMP {

// grammar for reading in files in the format of the Dual Decomposition algorithm of Torresani, Kolmogorov and Rother
namespace TorresaniEtAlInput {

   /* file format
   // Angular parentheses mean that it should be replaced with an integer number,
   // curly parentheses mean a floating point number.
   // Point and assignment id's are integers starting from 0.

   c comment line
   p <N0> <N1> <A> <E>     // # points in the left image, # points in the right image, # assignments, # edges
   a <a> <i0> <i1> {cost}  // specify assignment
   e <a> <b> {cost}        // specify edge

   i0 <id> {xi} {yi}       // optional - specify coordinate of a point in the left image
   i1 <id> {xi} {yi}       // optional - specify coordinate of a point in the left image
   n0 <i> <j>              // optional - specify that points <i> and <j> in the left image are neighbors
   n1 <i> <j>              // optional - specify that points <i> and <j> in the right image are neighbors
   */ 

   using parsing::mand_whitespace;
   using parsing::opt_whitespace;
   using parsing::positive_integer;
   using parsing::real_number;


   // first two integers are number of left nodes, number of right nodes, then comes number of assignments, and then number of quadratic terms
   struct no_left_nodes : pegtl::seq< positive_integer > {};
   struct no_right_nodes : pegtl::seq< positive_integer > {};
   struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'p'>, mand_whitespace, no_left_nodes, mand_whitespace, no_right_nodes, mand_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace > {};
   // numbers mean: assignment number (consecutive), then comes left node number, right node number, cost
   struct assignment : pegtl::seq < positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {};
   struct assignment_line : pegtl::seq< opt_whitespace, pegtl::string<'a'>, mand_whitespace, assignment, opt_whitespace> {};
   // numbers mean: number of left assignment, number of right assignment, cost
   struct quadratic_pot : pegtl::seq< positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number > {};
   struct quadratic_pot_line : pegtl::seq<opt_whitespace, pegtl::string<'e'>, mand_whitespace, quadratic_pot, opt_whitespace > {};

   struct comment_line : pegtl::seq< opt_whitespace, pegtl::string<'c'>, pegtl::until< pegtl::eol >> {};

   // artifacts from dual decomposition file format. We do not make use of them.
   struct neighbor_line : pegtl::seq< pegtl::sor<pegtl::string<'n','0'>, pegtl::string<'n','1'>>, pegtl::until< pegtl::eol>> {};
   struct coordinate_line : pegtl::seq< pegtl::sor<pegtl::string<'i','0'>, pegtl::string<'i','1'>>, pegtl::until< pegtl::eol>> {};

   // better way to cope with comment lines? On each line there may be a comment
   struct grammar : pegtl::must<
                    pegtl::star<comment_line>,
                    init_line,pegtl::eol,
                    pegtl::star< pegtl::sor<
                       pegtl::seq<quadratic_pot_line,pegtl::eol>,
                       pegtl::seq<assignment_line,pegtl::eol>,
                       comment_line,
                       neighbor_line,
                       coordinate_line,
                       pegtl::seq<opt_whitespace, pegtl::eol>, 
                       opt_whitespace
                    > >
                    //, pegtl::eof
                    > {};

   template< typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<> struct action< no_left_nodes > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_input& gmInput)
      {
         gmInput.no_left_nodes_ = std::stoul(in.string());
      }
   };
    
   template<> struct action< no_right_nodes > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_input& gmInput)
      {
         gmInput.no_right_nodes_ = std::stoul(in.string());
      }
   };
    
   template<> struct action< assignment > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_input& gmInput)
      {
         std::istringstream iss(in.string());
         std::size_t assignment_no; iss >> assignment_no;
         std::size_t left_node; iss >> left_node;
         std::size_t right_node; iss >> right_node;
         double cost; iss >> cost;

         assert(left_node < gmInput.no_left_nodes_);
         assert(right_node < gmInput.no_right_nodes_);
         assert(assignment_no == gmInput.assignment_.size());

         gmInput.assignment_.push_back({left_node,right_node,cost}); 
      }
   };
   template<> struct action< quadratic_pot > {
      template<typename INPUT>
      static void apply(const INPUT & in, graph_matching_input& gmInput)
      {
         std::istringstream iss(in.string());
         std::size_t assignment1; iss >> assignment1;
         std::size_t assignment2; iss >> assignment2;
         double cost; iss >> cost;

         assert(assignment1 < gmInput.assignment_.size());
         assert(assignment2 < gmInput.assignment_.size());

         gmInput.quadratic_.push_back({assignment1, assignment2, cost});
      }
   };

   graph_matching_input parse_file(const std::string& filename)
   {
      graph_matching_input gmInput;
      pegtl::file_parser problem(filename);
      std::cout << "parsing " << filename << "\n";

      const bool ret = problem.parse< grammar, action >( gmInput );
      std::sort(gmInput.assignment_.begin(), gmInput.assignment_.end());

      if(!ret) {
         throw std::runtime_error("could not read file " + filename);
      }

      return gmInput;
   }

   graph_matching_input parse_string(const std::string& input)
   {
      graph_matching_input gm_input;

      bool read_success = pegtl::parse<grammar, action>(input,"",gm_input);
      std::sort(gm_input.assignment_.begin(), gm_input.assignment_.end());

      if(!read_success) {
         throw std::runtime_error("could not read input");
      }

      return gm_input; 
   }

} // namespace TorresaniEtAlInput

} // namespace LPMP
