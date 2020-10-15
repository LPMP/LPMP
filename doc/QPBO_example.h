#pragma once
#include <vector>
#include "LP.h"
#include "vector.hxx"
#include "factors_messages.hxx"
#include "pegtl.hh"

namespace LPMP {

    using QPBO_instance = std::vector<std::vector<double>>;

    // parser
    //
    struct number : pegtl::seq< pegtl::opt< pegtl::string<'-'> >, pegtl::star< pegtl::blank >, pegtl::plus< pegtl::digit > > {};

    struct entry : pegtl::seq< pegtl::star< pegtl::blank >, number, pegtl::star< pegtl::blank >, pegtl::string<','> > {};
    struct row_begin : pegtl::seq<> {};
    struct row : pegtl::seq< row_begin, pegtl::star< entry >, pegtl::star< pegtl::blank >, number, pegtl::star< pegtl::blank >, pegtl::eolf > {};

    struct grammar : pegtl::plus< row > {};

    template< typename Rule >
        struct action
        : pegtl::nothing< Rule > {};

    template<> struct action< row_begin > { 
        template<typename Input>
            static void apply(const Input& in, QPBO_instance& input) 
            {   
                input.push_back( {} );
            }   
    };  

    template<> struct action< number > { 
        template<typename Input>
            static void apply(const Input& in, QPBO_instance& input) 
            {   
                input.back().push_back( std::stod(in.string()) );
            }   
    };  

    QPBO_instance parse_QPBO_file(const std::string filename)
    {
        QPBO_instance q;
        pegtl::file_parser problem(filename);

        const bool read_success = problem.parse<grammar, action>( q );
        return q; 
    }

    // factors

    class QPBO_diagonal_factor {
        public:
            QPBO_diagonal_factor(const double c) : cost(c) {}
            double LowerBound() const { return std::min(cost, 0.0); }
            double EvaluatePrimal() const { return primal*cost; }
            void init_primal() { primal = 0; }
            auto export_variables() { return std::tie(cost); }
            template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) {}
            template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) {}

            double cost;
            unsigned char primal;
    };

    class QPBO_off_diagonal_factor {
        public:
            QPBO_off_diagonal_factor(const double c) : cost({0.0, 0.0, c}) {}
            double LowerBound() const {
                return std::min({ 0.0, cost[0], cost[1], cost[0] + cost[1] + cost[2] });
            }
            double EvaluatePrimal() const {
                if( primal[0] * primal[1] != primal[2] )
                    return std::numeric_limits<double>::infinity();
                return primal[0]*cost[0] + primal[1]*cost[1] + primal[2]*cost[2];
            }
            void init_primal() { primal = {0,0,0}; }
            auto export_variables() { return std::tie(cost); }
            template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) {}
            template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) {}

            array<double,3> cost;
            std::array<unsigned char,3> primal;
    };

    // message

    class QPBO_message {
        public:
            QPBO_message(const unsigned char _index) : index(_index) {}
            void RepamLeft(QPBO_diagonal_factor& d, const double val, const std::size_t idx) { d.cost += val; }
            void RepamRight(QPBO_off_diagonal_factor& o, const double val, const std::size_t idx) { o.cost[index] += val; }
            template<typename MSG>
                void send_message_to_right(const QPBO_diagonal_factor& d, MSG& msg, const double scaling)
                {
                    msg[0] -= scaling*d.cost;
                }
            template<typename MSG>
                void send_message_to_left(const QPBO_off_diagonal_factor& o, MSG& msg, const double scaling)
                {
                    const double min_zero_marginal = std::min(0.0, o.cost[1-index]);
                    const double min_one_marginal = std::min(o.cost[index], o.cost[0] + o.cost[1] + o.cost[2]);
                    msg[0] -= scaling*(min_one_marginal - min_zero_marginal);
                }

        private:
            const unsigned char index;
    };

    // problem constructor

    template<typename FMC, typename DIAGONAL_FACTOR_CONTAINER, typename OFF_DIAGONAL_FACTOR_CONTAINER, typename QPBO_MESSAGE>
        class QPBO_problem_constructor {
            public:

                template<typename SOLVER>
                    QPBO_problem_constructor(SOLVER& s) : lp(&s.GetLP()) {}

                void construct (const QPBO_instance& q)
                {
                    const std::size_t n = q.size();
                    std::vector<DIAGONAL_FACTOR_CONTAINER*> diagonal_factors;
                    for(std::size_t i=0; i<n; ++i)
                        diagonal_factors.push_back( lp->template add_factor<DIAGONAL_FACTOR_CONTAINER>(q[i][i]) );
                    for(std::size_t i=0; i<n; ++i) {
                        for(std::size_t j=i+1; j<n; ++j) {
                            auto* o = lp->template add_factor<OFF_DIAGONAL_FACTOR_CONTAINER>(q[i][j] + q[j][i]);
                            lp->template add_message<QPBO_MESSAGE>(diagonal_factors[i], o, 0);
                            lp->template add_message<QPBO_MESSAGE>(diagonal_factors[j], o, 1);
                        } 
                    }
                }

            private:
                LP<FMC>* lp;
        };

    // FMC

    struct QPBO_FMC {
        constexpr static const char* name = "QPBO example"; 

        using diagonal_factor_container = FactorContainer<QPBO_diagonal_factor, QPBO_FMC, 0>;
        using off_diagonal_factor_container = FactorContainer<QPBO_off_diagonal_factor, QPBO_FMC, 1>;

        using QPBO_message_container = MessageContainer<QPBO_message, 0, 1, message_passing_schedule::left, variableMessageNumber, 2, QPBO_FMC, 0>;

        using FactorList = meta::list<diagonal_factor_container, off_diagonal_factor_container>;
        using MessageList = meta::list<QPBO_message_container>;

        using problem_constructor = QPBO_problem_constructor<QPBO_FMC, diagonal_factor_container, off_diagonal_factor_container, QPBO_message_container>;

    };

}

