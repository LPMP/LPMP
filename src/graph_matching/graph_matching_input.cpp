#include "graph_matching/graph_matching_input.h"
#include "pegtl_parse_rules.h"
#include <cassert>
#include "config.hxx"

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
   struct assignmentsline : pegtl::seq< opt_whitespace, pegtl::string<'a'>, mand_whitespace, assignment, opt_whitespace> {};
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
                       pegtl::seq<assignmentsline,pegtl::eol>,
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
         gmInput.no_left_nodes = std::stoul(in.string());
      }
   };
    
   template<> struct action< no_right_nodes > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_input& gmInput)
      {
         gmInput.no_right_nodes = std::stoul(in.string());
      }
   };
    
   template<> struct action< assignment > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_input& gmInput)
      {
         std::istringstream iss(in.string());
         std::size_t assignmentsno; iss >> assignmentsno;
         std::size_t left_node; iss >> left_node;
         std::size_t right_node; iss >> right_node;
         double cost; iss >> cost;

         if(left_node >= gmInput.no_left_nodes && left_node != graph_matching_input::no_assignment)
            throw std::runtime_error("error in assignment: left index larger than number of left nodes.");
         if(right_node >= gmInput.no_right_nodes && right_node != graph_matching_input::no_assignment)
            throw std::runtime_error("error in assignment: right index larger than number of right nodes.");
         if(assignmentsno != gmInput.assignments.size())
            throw std::runtime_error("error in assignment: assignment numbers must be sequential starting at 0."); 

         gmInput.add_assignment(left_node, right_node, cost);
      }
   };
   template<> struct action< quadratic_pot > {
      template<typename INPUT>
      static void apply(const INPUT & in, graph_matching_input& gmInput)
      {
         std::istringstream iss(in.string());
         std::size_t assignments1; iss >> assignments1;
         std::size_t assignments2; iss >> assignments2;
         double cost; iss >> cost;

         if(assignments1 >= gmInput.assignments.size())
            throw std::runtime_error("assignment number in quadratic infeasible: too large");
         if(assignments2 >= gmInput.assignments.size())
            throw std::runtime_error("assignment number in quadratic infeasible: too large");

         // check that quadratic assingment is feasible
         const std::size_t left_1 = gmInput.assignments[assignments1].left_node;
         const std::size_t right_1 = gmInput.assignments[assignments1].right_node;
         const std::size_t left_2 = gmInput.assignments[assignments2].left_node;
         const std::size_t right_2 = gmInput.assignments[assignments2].right_node;

         if(left_1 == left_2 && left_1 != graph_matching_input::no_assignment)
            throw std::runtime_error("assignments in quadratic infeasible: origin from same node");
         if(right_1 == right_2 && right_1 != graph_matching_input::no_assignment)
            throw std::runtime_error("assignments in quadratic infeasible: point to same node");

         gmInput.quadratic_terms.push_back({assignments1, assignments2, cost});
      }
   };

   graph_matching_input parse_file(const std::string& filename)
   {
      graph_matching_input gmInput;
      pegtl::file_parser problem(filename);
      if(debug())
          std::cout << "parsing " << filename << "\n";

      const bool ret = problem.parse< grammar, action >( gmInput );
      //std::sort(gmInput.assignments.begin(), gmInput.assignments.end());

      if(!ret) {
         throw std::runtime_error("could not read file " + filename);
      }

      return gmInput;
   }

   graph_matching_input parse_string(const std::string& input)
   {
      graph_matching_input gm_input;

      bool read_success = pegtl::parse<grammar, action>(input,"",gm_input);
      //std::sort(gm_input.assignments.begin(), gm_input.assignments.end());

      if(!read_success)
         throw std::runtime_error("could not read input");

      return gm_input; 
   }

} // namespace TorresaniEtAlInput

using parsing::mand_whitespace;
using parsing::opt_whitespace;
using parsing::positive_integer;

using partial_matching = std::vector<std::array<std::size_t,2>>;
struct assignment : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace, pegtl::string<'-','>'>, opt_whitespace, positive_integer, opt_whitespace > {};
struct grammar : pegtl::star<pegtl::sor<assignment,opt_whitespace,pegtl::eol>,pegtl::eof> {}; 

template< typename Rule >
struct action
: pegtl::nothing< Rule > {};

template<> struct action< assignment > {
   template<typename INPUT>
      static void apply(const INPUT& in, partial_matching& matchings)
      {
         std::istringstream iss(in.string());
         std::size_t left_node; iss >> left_node;
         char c; iss >> c;
         assert(c == '-');
         iss >> c;
         assert(c == '>');
         std::size_t right_node; iss >> left_node;

         matchings.push_back({left_node, right_node});
      }
};

graph_matching_input::labeling transform_partial_matching(partial_matching pm)
{
   const std::size_t no_nodes = (*std::max_element(pm.begin(), pm.end(), [](const auto& a, const auto& b) { return a[0] < b[0]; }))[0];
   graph_matching_input::labeling l;
   l.resize(no_nodes, std::numeric_limits<std::size_t>::max());
   for(const auto& m : pm)
      l[m[0]] = l[m[1]];
   return l;
}

graph_matching_input::labeling parse_graph_matching_result_file(const std::string& filename)
{
   partial_matching pm;
   throw std::runtime_error("not implemented yet");
   return transform_partial_matching(std::move(pm)); 
}

graph_matching_input::labeling parse_graph_matching_result_string(const std::string& input)
{
   partial_matching pm;
   const bool read_success = pegtl::parse<grammar, action>(input, "", pm);
   if(!read_success)
      throw std::runtime_error("could not read input");
   return transform_partial_matching(std::move(pm)); 
}

} // namespace LPMP
