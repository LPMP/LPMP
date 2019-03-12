#pragma once

#include "pegtl.hh"
#include "pegtl_parse_rules.h"
#include <cassert>
#include <string>

namespace LPMP {

    namespace cut_base_text_input {
 
        // import basic parsers
        using parsing::opt_whitespace;
        using parsing::mand_whitespace;
        using parsing::positive_integer;
        using parsing::real_number;

        struct empty_line : pegtl::seq<opt_whitespace, pegtl::eolf> {};
        struct comment_line : pegtl::seq< opt_whitespace, pegtl::sor<pegtl::string<'c'>, pegtl::string<'#'>>, pegtl::until<pegtl::eolf> > {};
        struct ignore_line : pegtl::sor<comment_line, empty_line> {};

        // #vertices #edges
        struct init_line : pegtl::opt< pegtl::seq< opt_whitespace, positive_integer, mand_whitespace, positive_integer, opt_whitespace, pegtl::eolf > > {};
        // $vertex_1 $vertex_2 $cost
        struct edge_line : pegtl::seq< opt_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace, pegtl::eolf > {};

        struct base_grammar : pegtl::seq< 
                              pegtl::star<ignore_line>,
                              init_line,
                              pegtl::star< pegtl::sor<ignore_line,edge_line> >
                              > {};

        template<typename IDENTIFIER>
            struct grammar : pegtl::seq< pegtl::star<ignore_line>, opt_whitespace, IDENTIFIER, opt_whitespace, pegtl::opt<pegtl::eol>, base_grammar > {};

        template< typename Rule >
            struct action
            : pegtl::nothing< Rule > {};

        template<> struct action< edge_line > {
            template<typename Input>
                static void apply(const Input& in, cut_base_instance& mc, const char sign) 
                {
                    std::istringstream iss(in.string());

                    std::size_t i;
                    iss >> i;

                    std::size_t j;
                    iss >> j;

                    double cost;
                    iss >> cost;

                    mc.add_edge(i, j, sign*cost);
                }
        };

        template<typename CUT_INSTANCE, typename IDENTIFIER, char SIGN = 1>
            CUT_INSTANCE parse_file(const std::string& filename)
            {
                CUT_INSTANCE input;
                pegtl::file_parser problem(filename);
                const bool read_success = problem.parse< grammar<IDENTIFIER>, action >(input, SIGN);
                if(!read_success) {
                    throw std::runtime_error("could not read file " + filename);
                }
                input.shift_to_zero_offset();

                return input;
            }

        template<typename CUT_INSTANCE, typename IDENTIFIER, char SIGN = 1>
            CUT_INSTANCE parse_string(const std::string& string)
            {
                CUT_INSTANCE input;
                const bool read_success = pegtl::parse<grammar<IDENTIFIER>, action >(string,"",input, SIGN);
                if(!read_success) {
                    throw std::runtime_error("could not read string");
                }
                input.shift_to_zero_offset();

                return input;
            }

    }

}
