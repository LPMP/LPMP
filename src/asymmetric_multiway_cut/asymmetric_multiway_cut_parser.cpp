#include "asymmetric_multiway_cut/asymmetric_multiway_cut_parser.h"
#include "pegtl_parse_rules.h"

namespace LPMP {

    namespace asymmetric_multiway_cut_parser {

        using parsing::opt_whitespace;
        using parsing::mand_whitespace;
        using parsing::positive_integer;
        using parsing::real_number;

struct empty_line : pegtl::seq<opt_whitespace, pegtl::eolf> {};
        struct comment_line : pegtl::seq< opt_whitespace, pegtl::sor<pegtl::string<'c'>, pegtl::string<'#'>>, pegtl::until<pegtl::eolf> > {};
        struct ignore_line : pegtl::sor<comment_line, empty_line> {};

        struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'A','S','Y','M','M','E','T','R','I','C',' ','M','U','L','T','I','W','A','Y',' ','C','U','T'>, opt_whitespace, pegtl::eolf> {};

        struct multicut_init_line : pegtl::seq< opt_whitespace, pegtl::string<'M','U','L','T','I','C','U','T'>, opt_whitespace, pegtl::eolf> {};
        struct edge_line : pegtl::seq< opt_whitespace, positive_integer, mand_whitespace, positive_integer, mand_whitespace, real_number, opt_whitespace, pegtl::eolf > {};

        struct node_costs_init_line : pegtl::seq< opt_whitespace, pegtl::string<'N','O','D','E',' ','C','O','S','T','S'>, opt_whitespace, pegtl::eolf> {};
        struct node_cost_line : pegtl::seq<opt_whitespace, real_number, pegtl::star<mand_whitespace, real_number>, opt_whitespace, pegtl::eolf> {};

        struct grammar : pegtl::seq<
                         init_line,
                         multicut_init_line,
                         pegtl::star<pegtl::sor<edge_line, ignore_line>>,
                         node_costs_init_line,
                         pegtl::star<pegtl::sor<node_cost_line, ignore_line>>,
                         pegtl::eof
                         > {};

        template< typename Rule >
            struct action
            : pegtl::nothing< Rule > {};

        template<> struct action< edge_line > {
            template<typename Input>
                static void apply(const Input& in, asymmetric_multiway_cut_instance& instance) 
                {
                    std::istringstream iss(in.string());

                    std::size_t i;
                    iss >> i;

                    std::size_t j;
                    iss >> j;

                    double capacity;
                    iss >> capacity;

                    instance.edge_costs.add_edge(i,j,capacity);
                }
        };

        template<> struct action< node_cost_line > {
            template<typename Input>
                static void apply(const Input& in, asymmetric_multiway_cut_instance& instance) 
                {
                    std::vector<double> costs;
                    std::istringstream iss(in.string());
                    double x;
                    while(iss >> x)
                        costs.push_back(x);

                    instance.node_costs.push_back(costs.begin(), costs.end()); 
                }
        };

        asymmetric_multiway_cut_instance parse_file(const std::string& filename)
        {
            asymmetric_multiway_cut_instance instance;
            pegtl::file_parser problem(filename);
            const bool read_success = problem.parse< grammar, action >(instance);
            if(!read_success) {
                throw std::runtime_error("could not read file " + filename);
            }

            return instance;
        }

        asymmetric_multiway_cut_instance parse_string(const std::string& string)
        {
            asymmetric_multiway_cut_instance instance;
            const bool read_success = pegtl::parse<grammar, action>(string,"",instance);
            if(!read_success) {
                throw std::runtime_error("could not read string");
            }

            return instance;
        }

    }

}
