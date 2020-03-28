#pragma once

#include <type_traits>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <fstream>
#include <sstream>

#include "LP.h"
#include "function_existence.hxx"
#include "template_utilities.hxx"
#include "tclap/CmdLine.h"

namespace LPMP {

    static std::vector<std::string> default_solver_options = {
        {""}, 
        {"-i"}, {""}, 
        {"--maxIter"}, {"1000"}
    };

    template<typename PROBLEM_CONSTRUCTOR>
        class primal_storage {
            public:
                template<typename T>
                    using has_get_primal_t = decltype(&T::get_primal);

                constexpr static bool can_get_primal();

                static auto get_primal_return_type()
                {
                    if constexpr(can_get_primal())
                        return static_cast<PROBLEM_CONSTRUCTOR*>(nullptr)->get_primal();
                    else
                        return nullptr;
                }

                void update_current_primal(PROBLEM_CONSTRUCTOR& pc);

                decltype(get_primal_return_type()) primal_solution_;
        };

    // class containing the LP, problem constructor list, input function and visitor
    // binds together problem constructors and solver and organizes input/output
    // base class for solvers with primal rounding, e.g. LP-based rounding heuristics, message passing rounding and rounding provided by problem constructors.
    template<typename LP_TYPE, typename VISITOR>
        class Solver : private primal_storage<typename LP_TYPE::FMC::problem_constructor> {

            public:
                using SolverType = Solver<LP_TYPE, VISITOR>;
                using FMC = typename LP_TYPE::FMC;
                using problem_constructor_type = typename FMC::problem_constructor;

                // default parameters

                Solver() : Solver(default_solver_options) {}

                Solver(int argc, char** argv);
                Solver(std::vector<std::string> options);

            private:
                Solver(const bool tmp);

            public:

                TCLAP::CmdLine& get_cmd() { return cmd_; }

                void Init_();

                const std::string& get_input_file() const { return inputFile_; }

                LPMP_FUNCTION_EXISTENCE_CLASS(has_order_factors, order_factors)
                    constexpr static bool can_order_factors();
                void order_factors();

                LPMP_FUNCTION_EXISTENCE_CLASS(HasWritePrimal,WritePrimal)
                    constexpr static bool CanWritePrimalIntoFile();
                constexpr static bool CanWritePrimalIntoString();

                const std::string& get_primal_string() const { return solution_; }

                void WritePrimal();

                std::string write_primal_into_string();

                LPMP_FUNCTION_EXISTENCE_CLASS(HasCheckPrimalConsistency,CheckPrimalConsistency)
                    constexpr static bool CanCheckPrimalConsistency();
                bool CheckPrimalConsistency() const;

                LPMP_FUNCTION_EXISTENCE_CLASS(HasTighten,Tighten)
                    constexpr static bool CanTighten();
                std::size_t Tighten(const std::size_t max_constraints_to_add);

                problem_constructor_type& GetProblemConstructor() { return problem_constructor_; }
                const problem_constructor_type& GetProblemConstructor() const { return problem_constructor_; }

                LP_TYPE& GetLP() { return lp_; }

                LPMP_FUNCTION_EXISTENCE_CLASS(has_solution,solution)
                    constexpr static bool visitor_has_solution();

                int Solve();

                // TODO: renamce Begin functino to pre_optimization or similar
                LPMP_FUNCTION_EXISTENCE_CLASS(has_begin,Begin)
                    constexpr static bool has_begin();

                // called before first iterations
                virtual void Begin();

                LPMP_FUNCTION_EXISTENCE_CLASS(has_pre_iterate, pre_iterate);
                constexpr static bool has_pre_iterate();

                // what to do before improving lower bound, e.g. setting reparametrization mode
                virtual void PreIterate(LpControl c);

                // what to do for improving lower bound, typically ComputePass or ComputePassAndPrimal
                virtual void Iterate(LpControl c);

                // what to do after one iteration of message passing, e.g. primal computation and/or tightening
                virtual void PostIterate(LpControl c);

                LPMP_FUNCTION_EXISTENCE_CLASS(HasEnd,End) 
                    // called after last iteration
                    constexpr static bool CanCallEnd();

                virtual void End();

                double EvaluatePrimal() const;
                // evaluate and register primal solution
                void RegisterPrimal();

                auto get_primal() const { return this->primal_solution_; }

                double lower_bound() const { return lowerBound_; }
                double primal_cost() const { return bestPrimalCost_; }

            protected:
                TCLAP::CmdLine cmd_;
                LP_TYPE lp_;
                problem_constructor_type problem_constructor_;
                primal_storage<problem_constructor_type> primal_storage_;

                // command line arguments
                TCLAP::ValueArg<std::string> inputFileArg_;
                TCLAP::ValueArg<std::string> outputFileArg_;
                std::string inputFile_;
                std::string outputFile_;
                TCLAP::ValueArg<INDEX> verbosity_arg_;

                double lowerBound_;
                // while Solver does not know how to compute primal, derived solvers do know. After computing a primal, they are expected to register their primals with the base solver
                double bestPrimalCost_ = std::numeric_limits<double>::infinity();
                std::string solution_;

                VISITOR visitor_;
                std::size_t iter = 0;
        };

    // local rounding interleaved with message passing 
    template<typename SOLVER>
        class MpRoundingSolver : public SOLVER
    {
        public:
            using SOLVER::SOLVER;

            virtual void Iterate(LpControl c);
    };

    // rounding based on primal heuristics provided by problem constructor
    template<typename SOLVER>
        class ProblemConstructorRoundingSolver : public SOLVER
    {
        public:
            using SOLVER::SOLVER;

            LPMP_FUNCTION_EXISTENCE_CLASS(HasComputePrimal,ComputePrimal)

                constexpr static bool can_compute_primal();

            void ComputePrimal();

            virtual void Begin();

            virtual void PostIterate(LpControl c);

            virtual void End();
    };

    // rounding based on (i) interleaved message passing followed by (ii) problem constructor rounding.
    // It is expected that (ii) takes (i)'s rounding into account.
    template<typename SOLVER>
        class CombinedMPProblemConstructorRoundingSolver : public ProblemConstructorRoundingSolver<SOLVER>
    {
        public:
            using ProblemConstructorRoundingSolver<SOLVER>::ProblemConstructorRoundingSolver;

            virtual void Iterate(LpControl c);

        private:
            std::size_t iter = 0;
            Direction cur_primal_computation_direction_ = Direction::forward; 
    };

    ////////////////////
    // implementation //
    ////////////////////

    template<typename PROBLEM_CONSTRUCTOR> 
        constexpr bool primal_storage<PROBLEM_CONSTRUCTOR>::can_get_primal()
        {
            return is_detected<has_get_primal_t, PROBLEM_CONSTRUCTOR>::value;
        }

    template<typename PROBLEM_CONSTRUCTOR> 
        void primal_storage<PROBLEM_CONSTRUCTOR>::update_current_primal(PROBLEM_CONSTRUCTOR& pc)
        {
            if constexpr(can_get_primal())
                primal_solution_ = pc.get_primal();
        }

    // class containing the LP, problem constructor list, input function and visitor
    // binds together problem constructors and solver and organizes input/output
    // base class for solvers with primal rounding, e.g. LP-based rounding heuristics, message passing rounding and rounding provided by problem constructors.

    template<typename LP_TYPE, typename VISITOR>
        Solver<LP_TYPE, VISITOR>::Solver(int argc, char** argv) : Solver(true) 
    {
        cmd_.parse(argc,argv);
        Init_(); 
    }

    template<typename LP_TYPE, typename VISITOR>
        Solver<LP_TYPE, VISITOR>::Solver(std::vector<std::string> options) : Solver(true)
    {
        cmd_.parse(options);
        Init_(); 
    }

    template<typename LP_TYPE, typename VISITOR>
        Solver<LP_TYPE, VISITOR>::Solver(const bool tmp)
        :
            cmd_(std::string("Command line options for ") + FMC::name, ' ', "0.0.1"),
            lp_(cmd_),
            inputFileArg_("i","inputFile","file from which to read problem instance",false,"","file name",cmd_),
            outputFileArg_("o","outputFile","file to write solution",false,"","file name",cmd_),
            verbosity_arg_("v","verbosity","verbosity level: 0 = silent, 1 = important runtime information, 2 = further diagnostics",false,1,"0,1,2",cmd_),
            visitor_(cmd_),
            problem_constructor_(*this)
    {
        std::cout << std::setprecision(10);
    }

    template<typename LP_TYPE, typename VISITOR>
        void Solver<LP_TYPE, VISITOR>::Init_()
        {
            try {  
                inputFile_ = inputFileArg_.getValue();
                outputFile_ = outputFileArg_.getValue();
                verbosity = verbosity_arg_.getValue();
                if(verbosity > 2) { throw TCLAP::ArgException("verbosity must be 0,1 or 2"); }
            } catch (TCLAP::ArgException &e) {
                std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
                exit(1);
            }
        }

    template<typename LP_TYPE, typename VISITOR>
        constexpr bool Solver<LP_TYPE, VISITOR>::can_order_factors()
        {
            return has_order_factors<problem_constructor_type, void>(); 
        }

    template<typename LP_TYPE, typename VISITOR>
        void Solver<LP_TYPE, VISITOR>::order_factors()
        {
            if constexpr(can_order_factors())
                problem_constructor_.order_factors();
        }

    template<typename LP_TYPE, typename VISITOR>
        constexpr bool Solver<LP_TYPE, VISITOR>::CanWritePrimalIntoFile()
        {
            return HasWritePrimal<problem_constructor_type, void, std::ofstream>();
        }

    template<typename LP_TYPE, typename VISITOR>
        constexpr bool Solver<LP_TYPE, VISITOR>::CanWritePrimalIntoString()
        {
            return HasWritePrimal<problem_constructor_type, void, std::stringstream>();
        } 

    template<typename LP_TYPE, typename VISITOR>
        void Solver<LP_TYPE, VISITOR>::WritePrimal()
        {
            if(outputFileArg_.isSet()) {
                std::ofstream output_file;
                output_file.open(outputFile_, std::ofstream::out);
                if(!output_file.is_open()) {
                    throw std::runtime_error("could not open file " + outputFile_);
                }

                output_file << solution_; // to do: not so nice, better have primal serialized, write back into factors, then write primal?

                //for_each_tuple(this->problemConstructor_, [this,&output_file](auto* l) {
                //      using pc_type = typename std::remove_pointer<decltype(l)>::type;
                //      static_if<SolverType::CanWritePrimalIntoFile<pc_type>()>([&](auto f) {
                //            f(*l).WritePrimal(output_file);
                //      });
                //}); 
            }
        }

    template<typename LP_TYPE, typename VISITOR>
        std::string Solver<LP_TYPE, VISITOR>::write_primal_into_string()
        {
            std::stringstream ss;

            if constexpr(CanWritePrimalIntoString())
                problem_constructor_.WritePrimal(ss);

            std::string sol = ss.str();
            return sol;
        }

    template<typename LP_TYPE, typename VISITOR>
        constexpr bool Solver<LP_TYPE, VISITOR>::CanCheckPrimalConsistency()
        {
            return HasCheckPrimalConsistency<problem_constructor_type, bool>();
        }

    template<typename LP_TYPE, typename VISITOR>
        bool Solver<LP_TYPE, VISITOR>::CheckPrimalConsistency() const
        {
            bool feasible = true;
            if constexpr(CanCheckPrimalConsistency())
                feasible = problem_constructor_.CheckPrimalConsistency();

            if(feasible)
                feasible = this->lp_.CheckPrimalConsistency();

            return feasible;
        }


    template<typename LP_TYPE, typename VISITOR>
        constexpr bool Solver<LP_TYPE, VISITOR>::CanTighten()
        {
            return HasTighten<problem_constructor_type, INDEX, INDEX>();
        }

    template<typename LP_TYPE, typename VISITOR>
        std::size_t Solver<LP_TYPE, VISITOR>::Tighten(const std::size_t max_constraints_to_add) 
        {
            INDEX constraints_added = 0;
            if constexpr(CanTighten())
                constraints_added += problem_constructor_.Tighten(max_constraints_to_add);
            return constraints_added;
        }

    template<typename LP_TYPE, typename VISITOR>
        constexpr bool Solver<LP_TYPE, VISITOR>::visitor_has_solution()
        {
            return has_solution<VISITOR, void, std::string>();
        } 

    template<typename LP_TYPE, typename VISITOR>
        int Solver<LP_TYPE, VISITOR>::Solve()
        {
            if(debug()) {
                std::cout << "lower bound before optimization = " << lp_.LowerBound() << "\n";
            }

            LpControl c = visitor_.begin(this->lp_);
            this->Begin();
            while(!c.end && !c.error) {
                this->PreIterate(c);
                this->Iterate(c);
                this->PostIterate(c);
                c = visitor_.visit(c, this->lowerBound_, this->bestPrimalCost_);
                ++iter;
            }
            if(!c.error) {
                this->End();
                RegisterPrimal();
                lowerBound_ = lp_.LowerBound();
                // possibly primal has been computed in end. Call visitor again
                visitor_.end(this->lowerBound_, this->bestPrimalCost_);
                if constexpr(visitor_has_solution()) {
                    this->visitor_.solution(this->solution_);
                }
                this->WritePrimal();
            }
            return !c.error;
        }


    template<typename LP_TYPE, typename VISITOR>
        constexpr bool Solver<LP_TYPE, VISITOR>::has_begin()
        {
            return has_begin<problem_constructor_type, void>();
        }

    template<typename LP_TYPE, typename VISITOR>
        void Solver<LP_TYPE, VISITOR>::Begin() 
        {
            lp_.Begin(); 
            if constexpr(has_begin())
                problem_constructor_.Begin();
            order_factors();
        }

    template<typename LP_TYPE, typename VISITOR>
        constexpr bool Solver<LP_TYPE, VISITOR>::has_pre_iterate()
        {
            return has_pre_iterate<problem_constructor_type, void>();
        }

    template<typename LP_TYPE, typename VISITOR>
        void Solver<LP_TYPE, VISITOR>::PreIterate(LpControl c) 
        {
            lp_.set_reparametrization(c.lp_repam);
            if constexpr(has_pre_iterate())
                problem_constructor_.pre_iterate();
        } 

    template<typename LP_TYPE, typename VISITOR>
        void Solver<LP_TYPE, VISITOR>::Iterate(LpControl c) 
        {
            lp_.ComputePass();
        } 

    template<typename LP_TYPE, typename VISITOR>
        void Solver<LP_TYPE, VISITOR>::PostIterate(LpControl c) 
        {
            if(c.computeLowerBound) {
                lowerBound_ = lp_.LowerBound();
                assert(std::isfinite(lowerBound_));
            }
            if(c.tighten) {
                Tighten(c.tightenConstraints);
            }
        } 

    template<typename LP_TYPE, typename VISITOR>
        constexpr bool Solver<LP_TYPE, VISITOR>::CanCallEnd()
        {
            return HasEnd<problem_constructor_type, void>();
        } 

    template<typename LP_TYPE, typename VISITOR>
        void Solver<LP_TYPE, VISITOR>::End() 
        {
            if constexpr(CanCallEnd())
                problem_constructor_.End();
            lp_.End();
        }

    template<typename LP_TYPE, typename VISITOR>
        double Solver<LP_TYPE, VISITOR>::EvaluatePrimal() const
        {
            if(CheckPrimalConsistency())
                return lp_.EvaluatePrimal();
            else
                return std::numeric_limits<double>::infinity();
        }

    template<typename LP_TYPE, typename VISITOR>
        void Solver<LP_TYPE, VISITOR>::RegisterPrimal()
        {
            const double cost = lp_.EvaluatePrimal();
            if(debug()) { std::cout << "register primal cost = " << cost << "\n"; }
            if(cost < bestPrimalCost_) {
                // assume solution is feasible
                const bool feasible = CheckPrimalConsistency();
                if(feasible) {
                    if(debug()) {
                        std::cout << "solution feasible\n";
                    }
                    bestPrimalCost_ = cost;
                    solution_ = write_primal_into_string();
                    this->update_current_primal(problem_constructor_);
                } else {
                    if(debug()) {
                        std::cout << "solution infeasible\n";
                    }
                } 
            }
        }


    template<typename SOLVER>
        void MpRoundingSolver<SOLVER>::Iterate(LpControl c)
        {
            if(c.computePrimal) {
                SOLVER::lp_.ComputeForwardPassAndPrimal();
                this->RegisterPrimal();
                SOLVER::lp_.ComputeBackwardPassAndPrimal();
                this->RegisterPrimal();
            } else {
                SOLVER::Iterate(c);
            }
        }

    template<typename SOLVER>
        constexpr bool ProblemConstructorRoundingSolver<SOLVER>::can_compute_primal()
        {
            return HasComputePrimal<typename SOLVER::problem_constructor_type, void>();
        }

    template<typename SOLVER>
        void ProblemConstructorRoundingSolver<SOLVER>::ComputePrimal()
        {
            SOLVER::lp_.init_primal();
            // compute the primal in parallel.
            // for this, first we have to wait until the rounding procedure has read off everything from the LP model before optimizing further
            if constexpr(can_compute_primal())
                this->problem_constructor_.ComputePrimal();
        }

    template<typename SOLVER>
        void ProblemConstructorRoundingSolver<SOLVER>::Begin()
        {
            SOLVER::Begin();
            //ComputePrimal();
            //this->RegisterPrimal();
        }

    template<typename SOLVER>
        void ProblemConstructorRoundingSolver<SOLVER>::PostIterate(LpControl c)
        {
            if(c.computePrimal) {
                ComputePrimal();
                this->RegisterPrimal();
            }
            SOLVER::PostIterate(c);
        }

    template<typename SOLVER>
        void ProblemConstructorRoundingSolver<SOLVER>::End()
        {
            SOLVER::End(); // first let problem constructors end (done in Solver)
            this->RegisterPrimal();
        }

    template<typename SOLVER>
        void CombinedMPProblemConstructorRoundingSolver<SOLVER>::Iterate(LpControl c)
        {
            if(c.computePrimal) {
                this->RegisterPrimal();
                // alternatively  compute forward and backward based rounding
                if(cur_primal_computation_direction_ == Direction::forward) {
                    if(verbosity >= 2) { std::cout << "compute primal for forward pass\n"; }
                    this->lp_.ComputeForwardPassAndPrimal();
                    this->lp_.ComputeBackwardPass();
                    cur_primal_computation_direction_ = Direction::backward;
                } else {
                    assert(cur_primal_computation_direction_ == Direction::backward);
                    if(verbosity >= 2) { std::cout << "compute primal for backward pass\n"; }
                    this->lp_.ComputeForwardPass();
                    this->lp_.ComputeBackwardPassAndPrimal();
                    cur_primal_computation_direction_ = Direction::forward;
                }
            } else {
                SOLVER::Iterate(c);
            }
            ++iter;
        }

} // end namespace LPMP
