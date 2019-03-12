#include "multicut/multicut_text_input.h"
#include "pegtl.hh"
#include "pegtl_parse_rules.h"
#include <cassert>

namespace LPMP {

namespace multicut_text_input {

// import basic parsers
using parsing::opt_whitespace;
using parsing::mand_whitespace;
using parsing::positive_integer;
using parsing::real_number;

struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'M','U','L','T','I','C','U','T'>, opt_whitespace, pegtl::eolf> {};
struct empty_line : pegtl::seq<opt_whitespace, pegtl::eolf> {};
struct comment_line : pegtl::seq< opt_whitespace, pegtl::sor<pegtl::string<'c'>, pegtl::string<'#'>>, pegtl::until<pegtl::eolf> > {};
struct ignore_line : pegtl::sor<comment_line, empty_line> {};

struct edge_line : pegtl::seq< opt_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace, pegtl::eolf > {};

struct grammar : pegtl::seq< init_line, pegtl::star< pegtl::sor<ignore_line,edge_line> > > {};

template< typename Rule >
struct action
: pegtl::nothing< Rule > {};

template<> struct action< edge_line > {
    template<typename Input>
        static void apply(const Input& in, multicut_instance& mc) 
        {
            std::istringstream iss(in.string());

            std::size_t i;
            iss >> i;

            std::size_t j;
            iss >> j;

            double capacity;
            iss >> capacity;

            mc.add_edge(i,j,capacity);
        }
};

multicut_instance parse_file(const std::string& filename)
{
       multicut_instance input;
       pegtl::file_parser problem(filename);
       const bool read_success = problem.parse< grammar, action >(input);
       if(!read_success) {
           throw std::runtime_error("could not read file " + filename);
       }

       return input;
}

multicut_instance parse_string(const std::string& string)
{
       multicut_instance input;
       const bool read_success = pegtl::parse<grammar, action>(string,"",input);
       if(!read_success) {
           throw std::runtime_error("could not read string");
       }

       return input;
}

} // namespace multicut_text_input



/*
namespace MulticutTextInput {

   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::positive_integer;
   using Parsing::real_number;

   // add pegtl::eolf to line structs and remove from grammar
   struct comment : pegtl::seq< opt_whitespace, pegtl::string<'#'>, pegtl::until< pegtl::eolf > > {};
   struct comment_line : pegtl::seq< pegtl::sor< comment, opt_whitespace >, pegtl::eolf > {};
   struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'M','U','L','T','I','C','U','T'>, opt_whitespace > {};
   struct numberOfVariables_line : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace > {};
   struct edge_line : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace, positive_integer, opt_whitespace, real_number, opt_whitespace, pegtl::eolf > {};

   struct grammar : pegtl::must<
                    init_line, pegtl::eol,
                    numberOfVariables_line, pegtl::eol,
                    pegtl::star<edge_line, pegtl::eol>,
                    pegtl::opt<edge_line>,
                    pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>,
                    pegtl::eof> {};

   struct lifted_line : pegtl::seq<opt_whitespace, pegtl::string<'L','I','F','T','E','D'>, opt_whitespace> {};
   struct lifted_edge_line : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace, positive_integer, opt_whitespace, real_number, opt_whitespace > {};

   struct LiftedMulticutGrammar : pegtl::must<
                    init_line, pegtl::eol,
                    numberOfVariables_line, pegtl::eol,
                    pegtl::star<edge_line, pegtl::eol>,
                    pegtl::opt<edge_line>,
                    pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>,
                    // now come lifted edges //
                    lifted_line, pegtl::eol,
                    pegtl::star<lifted_edge_line, pegtl::eol>,
                    pegtl::opt<lifted_edge_line>,
                    pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>,
                    pegtl::eof> {};

   struct base_edges_begin_line : pegtl::seq< opt_whitespace, pegtl::string<'B','A','S','E',' ','E','D','G','E','S'>, opt_whitespace, pegtl::eol > {};
   struct triplet_begin_line : pegtl::seq< opt_whitespace, pegtl::string<'T','R','I','P','L','E','T','S'>, opt_whitespace, pegtl::eol > {};
   struct triplet_line : pegtl::seq<
                         opt_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace,
                         real_number, mand_whitespace, real_number, mand_whitespace, real_number, mand_whitespace, real_number, mand_whitespace, real_number, opt_whitespace, 
                         pegtl::eolf> {};

   struct higher_order_grammar : pegtl::must<
                                 pegtl::star<comment_line>,
                                 base_edges_begin_line,
                                 pegtl::star< pegtl::sor< comment_line, edge_line > >,
                                 triplet_begin_line,
                                 pegtl::star< pegtl::sor< comment_line, triplet_line > >,
                                 pegtl::at< pegtl::eof >
                                 > {};


   template<typename SOLVER, typename Rule >
      struct action
      : pegtl::nothing< Rule > {};


   template<typename SOLVER> struct action<SOLVER, edge_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, SOLVER& mc)
      {
         std::stringstream s(in.string());
         INDEX i1; s >> i1;
         INDEX i2; s >> i2;
         REAL cost; s >> cost;
         if(i1 > i2) {
            std::swap(i1,i2);
         }
         assert(i1 < i2);
         mc.AddUnaryFactor( i1,i2,cost );
      }
   };

   template<typename SOLVER> struct action<SOLVER, triplet_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, SOLVER& mc)
      {
         std::stringstream s(in.string());
         INDEX u; s >> u;
         INDEX v; s >> v;
         INDEX w; s >> w;
         REAL cost_000; s >> cost_000;
         REAL cost_011; s >> cost_011;
         REAL cost_101; s >> cost_101;
         REAL cost_110; s >> cost_110;
         REAL cost_111; s >> cost_111;
         assert(u < v && v << w);
         mc.add_higher_order_triplet( u,v,w, cost_000, cost_011, cost_101, cost_110, cost_111 );
      }
   };


   template<typename SOLVER> struct action<SOLVER, lifted_edge_line > {
      template<typename INPUT>
      static void apply(const INPUT & in, SOLVER& mc)
      {
         std::stringstream s(in.string());
         INDEX i1; s >> i1;
         INDEX i2; s >> i2;
         REAL cost; s >> cost;
         if(i1 > i2) {
            std::swap(i1,i2);
         }
         assert(i1 < i2);
         mc.AddLiftedUnaryFactor( i1,i2,cost );
      }
   };

   template<typename SOLVER>
      struct actionSpecialization {
         template<typename RULE> struct type : public action<SOLVER,RULE> {};
      };

      
   template<typename SOLVER>
   bool ParseProblem(const std::string filename, SOLVER& pd)
   {
      std::cout << "parsing " << filename << "\n";

      pegtl::file_parser problem(filename);
      auto& mc = pd.template GetProblemConstructor<0>();
      return problem.parse< grammar, actionSpecialization<decltype(mc)>::template type >(mc);
   }

   template<typename SOLVER>
   bool ParseLiftedProblem(const std::string filename, SOLVER& pd)
   {
      std::cout << "parsing " << filename << "\n";

      pegtl::file_parser problem(filename);
      auto& mc = pd.template GetProblemConstructor<0>();
      return problem.parse< LiftedMulticutGrammar, actionSpecialization<decltype(mc)>::template type >(mc);
   }

   template<typename SOLVER>
   bool parse_higher_order(const std::string filename, SOLVER& s)
   {
      std::cout << "parsing " << filename << "\n";
      pegtl::file_parser problem(filename);
      auto& mc = s.template GetProblemConstructor<0>();
      auto ret = problem.parse< higher_order_grammar, actionSpecialization<decltype(mc)>::template type >(mc);
      if(verbosity >= 1) {
         std::cout << "loaded problem with " << mc.number_of_edges() << " edges and " << mc.number_of_triplets() << " triplets\n";
      }
      return ret;
   }


} // end namespace MulticutTextInput
*/


} // namespace LPMP
