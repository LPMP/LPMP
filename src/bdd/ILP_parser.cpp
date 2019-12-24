#include "bdd/ILP_parser.h"
#include "pegtl.hh"
#include "pegtl_parse_rules.h"
#include "bdd/ILP_input.h"

namespace LPMP {

    namespace ILP_parser {

        // import basic parsers
        using parsing::opt_whitespace;
        using parsing::mand_whitespace;
        using parsing::opt_invisible;
        using parsing::mand_invisible;
        using parsing::positive_integer;
        using parsing::real_number;

        struct min_line : pegtl::seq<opt_whitespace, pegtl::string<'M','i','n','i','m','i','z','e'>, opt_whitespace, pegtl::eol> {};

        struct sign : pegtl::sor<pegtl::string<'+'>, pegtl::string<'-'>> {};

        struct variable_name : pegtl::seq< 
                               pegtl::alpha, pegtl::star< pegtl::sor< pegtl::alnum, pegtl::string<'_'>, pegtl::string<'-'>, pegtl::string<'/'>, pegtl::string<'('>, pegtl::string<')'>, pegtl::string<'{'>, pegtl::string<'}'>, pegtl::string<','> > > 
                               > {};

        struct variable_coefficient : pegtl::seq< pegtl::opt<real_number> > {};

        struct variable : pegtl::seq<variable_coefficient, opt_whitespace, variable_name> {};

        struct weighted_sum_of_variables : pegtl::seq< opt_whitespace, variable, pegtl::star< opt_whitespace, sign, opt_whitespace, variable >>{};

        struct objective_coefficient : real_number {};
        struct objective_variable : variable_name {};
        struct objective_term : pegtl::seq< pegtl::opt<sign, opt_whitespace>, pegtl::opt<objective_coefficient, opt_whitespace, pegtl::opt<pegtl::string<'*'>>, opt_whitespace>, objective_variable> {};
        struct subject_to : pegtl::string<'S','u','b','j','e','c','t',' ','T','o'> {};
        struct objective_line : pegtl::seq< pegtl::not_at<subject_to>, pegtl::star<opt_whitespace, objective_term>, opt_whitespace, pegtl::eol> {};

        struct subject_to_line : pegtl::seq<opt_whitespace, subject_to, opt_whitespace, pegtl::eol> {};

        struct inequality_type : pegtl::sor<pegtl::string<'<','='>, pegtl::string<'>','='>, pegtl::string<'='>> {};

        struct inequality_identifier : pegtl::seq<pegtl::alpha, pegtl::star<pegtl::alnum>, opt_whitespace, pegtl::string<':'>, opt_whitespace> {};
        struct new_inequality : pegtl::seq<opt_whitespace, pegtl::not_at<pegtl::sor<pegtl::string<'E','n','d'>,pegtl::string<'B','o','u','n','d','s'>>>, pegtl::opt<inequality_identifier>> {};

        struct inequality_coefficient : pegtl::seq<pegtl::digit, pegtl::star<pegtl::digit>> {};
        struct inequality_variable : variable_name {};
        struct inequality_term : pegtl::seq< pegtl::opt<sign, opt_whitespace>, pegtl::opt<inequality_coefficient, opt_whitespace, pegtl::opt<pegtl::string<'*'>>, opt_whitespace>, inequality_variable> {};
        struct right_hand_side : pegtl::seq< pegtl::opt<sign>, pegtl::digit, pegtl::star<pegtl::digit> > {};

        struct inequality_line : pegtl::seq< new_inequality, 
            pegtl::star<opt_whitespace, inequality_term, opt_whitespace, pegtl::opt<pegtl::eol>>,
            opt_whitespace, inequality_type, opt_whitespace, right_hand_side, opt_whitespace, pegtl::eol> {};

        struct bounds_begin : pegtl::seq<opt_whitespace, pegtl::string<'B','o','u','n','d','s'>, opt_whitespace, pegtl::eol> {};
        struct binaries : pegtl::seq<opt_whitespace, pegtl::string<'B','i','n','a','r','i','e','s'>, opt_whitespace, pegtl::eol> {};
        struct binary_variable : variable_name {};
        struct list_of_binary_variables : pegtl::seq< 
            pegtl::star<opt_whitespace, pegtl::opt<pegtl::eol>, opt_whitespace, pegtl::not_at< pegtl::string<'E','n','d'> >, binary_variable>,
            opt_whitespace, pegtl::eol > {};

        struct bounds : pegtl::seq<pegtl::opt<bounds_begin, binaries, list_of_binary_variables> > {};

        struct end_line : pegtl::seq<opt_whitespace, pegtl::string<'E','n','d'>, opt_whitespace, pegtl::eolf> {};

        struct grammar : pegtl::seq<
                         min_line,
                         pegtl::star<objective_line>,
                         subject_to_line,
                         pegtl::star<inequality_line>,
                         bounds,
                         end_line> {};

        template< typename Rule >
            struct action
            : pegtl::nothing< Rule > {};

        struct tmp_storage {
            double objective_coeff = 1.0;
            int constraint_coeff = 1;
        };

        template<> struct action< sign > {
            template<typename INPUT>
                static void apply(const INPUT & in, ILP_input& i, tmp_storage& tmp)
                {
                    if(in.string() == "+") {
                    } else if(in.string() == "-") {
                        tmp.objective_coeff *= -1.0;
                        tmp.constraint_coeff *= -1;
                    } else
                        throw std::runtime_error("sign not recognized");
                }
        };

        template<> struct action< objective_coefficient > {
            template<typename INPUT>
                static void apply(const INPUT & in, ILP_input& i, tmp_storage& tmp)
                {
                    tmp.objective_coeff *= std::stod(in.string());
                }
        };

        template<> struct action< objective_variable > {
            template<typename INPUT>
                static void apply(const INPUT & in, ILP_input& i, tmp_storage& tmp)
                {
                    const std::string var = in.string();
                    i.add_to_objective(tmp.objective_coeff, var);
                    tmp = tmp_storage{};
                }
        };

        template<> struct action< new_inequality > {
            template<typename INPUT>
                static void apply(const INPUT & in, ILP_input& i, tmp_storage& tmp)
                {
                    i.begin_new_inequality();
                }
        };

        template<> struct action< inequality_coefficient > {
            template<typename INPUT>
                static void apply(const INPUT & in, ILP_input& i, tmp_storage& tmp)
                {
                    tmp.constraint_coeff *= std::stoi(in.string());
                }
        };

        template<> struct action< inequality_variable > {
            template<typename INPUT>
                static void apply(const INPUT & in, ILP_input& i, tmp_storage& tmp)
                {
                    const std::string var = in.string();
                    i.add_to_constraint(tmp.constraint_coeff, var);
                    tmp = tmp_storage{};
                }
        };

        template<> struct action< inequality_type > {
            template<typename INPUT>
                static void apply(const INPUT & in, ILP_input& i, tmp_storage& tmp)
                {
                    const std::string ineq = in.string();
                    if(ineq == "=")
                        i.set_inequality_type(LPMP::inequality_type::equal);
                    else if(ineq == "<=")
                        i.set_inequality_type(LPMP::inequality_type::smaller_equal);
                    else if(ineq == ">=")
                        i.set_inequality_type(LPMP::inequality_type::greater_equal);
                    else
                        throw std::runtime_error("inequality type not recognized");
                }
        };

        template<> struct action< right_hand_side > {
            template<typename INPUT>
                static void apply(const INPUT & in, ILP_input& i, tmp_storage& tmp)
                {
                    const double val = std::stod(in.string());
                    i.set_right_hand_side(val);
                }
        };

        ILP_input parse_file(const std::string& filename)
        {
            ILP_input input;
            tmp_storage tmp;
            pegtl::file_parser problem(filename);
            const bool success = problem.parse< grammar, action >(input, tmp); 
            if(!success)
                throw std::runtime_error("could not read input file " + filename);
            return input;
        }

        ILP_input parse_string(const std::string& input_string)
        {
            ILP_input input;
            tmp_storage tmp;

            const bool success = pegtl::parse<grammar, action>(input_string,"", input, tmp);
            if(!success)
                throw std::runtime_error("could not read input");
            return input;
        }

    } 

}
