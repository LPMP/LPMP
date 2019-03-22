#include "graph_matching/graph_matching_koopmans_beckmann_input.h"
#include "pegtl_parse_rules.h"
#include <cassert>

namespace LPMP {

namespace graph_matching_koopmans_beckmann_input {

   /* file format
   c comment line
   L
   x_11, ..., x_n1
   .
   .
   .
   x_1n, ..., x_nn

   A
   x_11, ..., x_n1
   .
   .
   .
   x_1n, ..., x_nn

   B
   x_11, ..., x_n1
   .
   .
   .
   x_1n, ..., x_nn
   */ 

   using parsing::mand_whitespace;
   using parsing::opt_whitespace;
   using parsing::positive_integer;
   using parsing::real_number;

   // read matrix in csv format
   struct matrix_input {
       std::vector<double> entries;
       std::size_t no_cols = std::numeric_limits<std::size_t>::max();

       void reset()
       {
           entries.clear();
           no_cols = std::numeric_limits<std::size_t>::max();
       }

       Eigen::MatrixXd matrix() const {
           if(entries.size() % no_cols != 0)
               throw std::runtime_error("dimensions of matrix not correct.");
           const std::size_t no_rows = entries.size() / no_cols;
           Eigen::MatrixXd m(no_cols, no_rows);
           std::size_t c=0;
           for(std::size_t i=0; i<m.cols(); ++i) {
               for(std::size_t j=0; j<m.rows(); ++j) {
                   m(i,j) = entries[c++];
               }
           }
           return m; 
       }
   };

   struct matrix_begin : pegtl::string<> {};
   struct matrix_entry : pegtl::seq< real_number > {};
   struct matrix_separator : pegtl::opt<pegtl::string<','>> {};
   struct matrix_row_end : pegtl::seq< pegtl::eolf > {};
   struct matrix : pegtl::seq< matrix_begin, 
    pegtl::star< opt_whitespace, matrix_entry, pegtl::until< matrix_row_end, pegtl::seq< opt_whitespace, matrix_separator, opt_whitespace, matrix_entry >>
                   > > {};

   struct comment_line : pegtl::seq< opt_whitespace, pegtl::string<'c'>, pegtl::until< pegtl::eolf >> {};
   struct empty_line : pegtl::seq< opt_whitespace, pegtl::eolf >{};
   struct ignore_line : pegtl::sor<comment_line, empty_line> {};

   struct L_end : pegtl::string<> {};
   struct L : pegtl::seq<pegtl::string<'L'>, opt_whitespace, pegtl::eol, matrix, L_end> {};

   struct A_end : pegtl::string<> {};
   struct A : pegtl::seq<pegtl::string<'A'>, opt_whitespace, pegtl::eol, matrix, A_end> {};

   struct B_end : pegtl::string<> {};
   struct B : pegtl::seq<pegtl::string<'B'>, opt_whitespace, pegtl::eol, matrix, B_end> {};

   struct qaplib_A : pegtl::seq<opt_whitespace, matrix, A_end> {};
   struct qaplib_B : pegtl::seq<opt_whitespace, matrix, B_end> {};

   // better way to cope with comment lines? On each line there may be a comment
   struct LPMP_grammar : pegtl::seq<
                    pegtl::star<ignore_line>,
                    L,
                    pegtl::star<ignore_line>,
                    A,
                    pegtl::star<ignore_line>,
                    B,
                    pegtl::star<ignore_line>
                    > {};

   struct qaplib_grammar : pegtl::seq<
                           pegtl::star<ignore_line>,
                           opt_whitespace, positive_integer, opt_whitespace, pegtl::eol,
                           pegtl::star<ignore_line>,
                           qaplib_A,
                           pegtl::star<ignore_line>,
                           qaplib_B,
                           pegtl::star<ignore_line>
                           > {};

   struct grammar : pegtl::sor<LPMP_grammar, qaplib_grammar> {};

   template< typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<> struct action< L_end > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_instance_koopmans_beckmann& gm, matrix_input& mi)
      {
         gm.L = mi.matrix();
      }
   };

   template<> struct action< A_end > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_instance_koopmans_beckmann& gm, matrix_input& mi)
      {
         gm.A = mi.matrix();
      }
   };

   template<> struct action< B_end > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_instance_koopmans_beckmann& gm, matrix_input& mi)
      {
         gm.B = mi.matrix();
      }
   };

   // linear cost not specified in qaplib format
   template<> struct action< qaplib_B > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_instance_koopmans_beckmann& gm, matrix_input& mi)
      {
         gm.L = Eigen::MatrixXd(gm.A.rows(), gm.A.cols());
         gm.L.setZero();
      }
   };

   template<> struct action< matrix_begin > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_instance_koopmans_beckmann& gm, matrix_input& mi)
      {
         mi.reset();
      }
   };

   template<> struct action< matrix_entry > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_instance_koopmans_beckmann& gm, matrix_input& mi)
      {
         const double x = std::stod(in.string());
         mi.entries.push_back(x);
      }
   };

   template<> struct action< matrix_row_end > {
      template<typename INPUT>
      static void apply(const INPUT& in, graph_matching_instance_koopmans_beckmann& gm, matrix_input& mi)
      {
          if(mi.no_cols == std::numeric_limits<std::size_t>::max())
              mi.no_cols = mi.entries.size();
          else if(mi.entries.size() % mi.no_cols != 0)
              throw std::runtime_error("dimensions of matrix not correct.");
      }
   };

   graph_matching_instance_koopmans_beckmann parse_file(const std::string& filename)
   {
      graph_matching_instance_koopmans_beckmann gm;
      matrix_input mi;
      pegtl::file_parser problem(filename);
      std::cout << "parsing " << filename << "\n";

      const bool ret = problem.parse< grammar, action >( gm, mi );
      //std::sort(gm.assignments.begin(), gm.assignments.end());

      if(!ret) {
         throw std::runtime_error("could not read file " + filename);
      }

      return gm;
   }

   graph_matching_instance_koopmans_beckmann parse_string(const std::string& input)
   {
      graph_matching_instance_koopmans_beckmann gm;
      matrix_input mi;

      bool read_success = pegtl::parse<grammar, action>(input,"",gm,mi);
      //std::sort(gm_input.assignments.begin(), gm_input.assignments.end());

      if(!read_success)
         throw std::runtime_error("could not read input");

      return gm; 
   }

}

} // namespace LPMP
