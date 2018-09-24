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

} // namespace LPMP
