#include "test_mrf.hxx"

using namespace LPMP;

int main()
{
    using FMC = FMC_SRMP_T;
    using VisitorType = StandardTighteningVisitor;
    using SolverType = MpRoundingSolver<Solver<LP<FMC>,VisitorType>>;

    auto tightening_solver_options = solver_options;
    tightening_solver_options.push_back("--tighten");
    tightening_solver_options.push_back("--tightenIteration");
    tightening_solver_options.push_back("5");
    tightening_solver_options.push_back("--tightenConstraintsMax");
    tightening_solver_options.push_back("1");
    tightening_solver_options.push_back("--tightenInterval");
    tightening_solver_options.push_back("10");
    tightening_solver_options.push_back("--tightenReparametrization");
    tightening_solver_options.push_back("uniform:0.5");


    // binary triplet
    {
        SolverType s(tightening_solver_options);
        auto& mrf = s.template GetProblemConstructor<0>();
        mrf.add_unary_factor(std::vector<REAL>(2,0.0));
        mrf.add_unary_factor(std::vector<REAL>(2,0.0));
        mrf.add_unary_factor(std::vector<REAL>(2,0.0));

        // make a cycle of length 3 visiting each label once -> one negative Potts and two positive Potts
        // make number of labels varying to check whether bounds work alright
        auto pos_potts = construct_potts(2,2, 0.0, 1.0);
        auto neg_potts = construct_potts(2,2, 1.0, 0.0);
        mrf.add_pairwise_factor(0,1,neg_potts);
        mrf.add_pairwise_factor(0,2,pos_potts);
        mrf.add_pairwise_factor(1,2,pos_potts);

        mrf.add_tightening_triplet(0,1,2);
        s.Solve();
        test(std::abs(s.lower_bound() - 1.0) <= eps);
    }

    // multi-label triplet
    {
        SolverType s(tightening_solver_options);
        auto& mrf = s.template GetProblemConstructor<0>();
        mrf.add_unary_factor(std::vector<REAL>(2,0.0));
        mrf.add_unary_factor(std::vector<REAL>(3,0.0));
        mrf.add_unary_factor(std::vector<REAL>(4,0.0));

        // make a cycle of length 3 visiting each label once -> one negative Potts and two positive Potts
        // make number of labels varying to check whether bounds work alright
        auto negPotts23 = construct_potts(2,3, 1.0, 0.0);
        auto posPotts24 = construct_potts(2,4, 0.0, 1.0);
        auto posPotts34 = construct_potts(3,4, 0.0, 1.0);
        mrf.add_pairwise_factor(0,1,negPotts23);
        mrf.add_pairwise_factor(0,2,posPotts24);
        mrf.add_pairwise_factor(1,2,posPotts34);

        mrf.add_tightening_triplet(0,1,2);
        s.Solve();
        test(std::abs(s.lower_bound() - 1.0) <= eps);
    }

    // triplet search
    {
        SolverType s(tightening_solver_options);
        auto& mrf = s.template GetProblemConstructor<0>(); 
        mrf.add_unary_factor(std::vector<REAL>(2,0.0));
        mrf.add_unary_factor(std::vector<REAL>(3,0.0));
        mrf.add_unary_factor(std::vector<REAL>(4,0.0));

        // make a cycle of length 3 visiting each label once -> one negative Potts and two positive Potts
        // make number of labels varying to check whether bounds work alright
        auto neg_potts_23 = construct_potts(2,3, 1.0, 0.0);
        auto pos_potts_24 = construct_potts(2,4, 0.0, 1.0);
        auto pos_potts_34 = construct_potts(3,4, 0.0, 1.0);
        mrf.add_pairwise_factor(0,1,neg_potts_23);
        mrf.add_pairwise_factor(0,2,pos_potts_24);
        mrf.add_pairwise_factor(1,2,pos_potts_34);

        s.Solve();
        test(std::abs(s.lower_bound() - 1.0) <= eps);
    }

    // cycle: k projection graph (individual labels) can find violated cycle
    {
        SolverType s(tightening_solver_options);
        auto& mrf = s.template GetProblemConstructor<0>();
        mrf.add_unary_factor(std::vector<REAL>(2,0.0));
        mrf.add_unary_factor(std::vector<REAL>(2,0.0));
        mrf.add_unary_factor(std::vector<REAL>(2,0.0));
        mrf.add_unary_factor(std::vector<REAL>(2,0.0));

        // make a cycle of length 4 visiting each label once -> one negative Potts and three positive Potts
        auto pos_potts = construct_potts(2,2, 0.0, 1.0);
        auto neg_potts = construct_potts(2,2, 1.0, 0.0);
        mrf.add_pairwise_factor(0,1,neg_potts);
        mrf.add_pairwise_factor(1,2,pos_potts);
        mrf.add_pairwise_factor(2,3,pos_potts);
        mrf.add_pairwise_factor(0,3,pos_potts);

        s.Solve();
        test(std::abs(s.lower_bound() - 1.0) <= eps);
    }

    // to do: add test problem where only create_expanded_projection_graph, but not create_k_projection_graph, can find a violated cycle.
    // to do: modified cycle inequalites do not work -> correct
    /*
    {
       // cycle: expanded projection graph
       SolverType s(tightening_solver_options);
       auto& mrf = s.template GetProblemConstructor<0>();
       mrf.add_unary_factor(std::vector<REAL>(4,0.0));
       mrf.add_unary_factor(std::vector<REAL>(4,0.0));
       mrf.add_unary_factor(std::vector<REAL>(4,0.0));
       mrf.add_unary_factor(std::vector<REAL>(4,0.0));

       std::vector<REAL> Potts2 = {0,0,1,1, 0,0,1,1, 1,1,0,0, 1,1,0,0};
       std::vector<REAL> negPotts2;
       std::transform(Potts2.begin(), Potts2.end(), std::back_inserter(negPotts2), [](auto x) { return -x+1; });

       mrf.add_pairwise_factor(0,1,negPotts2);
       mrf.add_pairwise_factor(1,2,Potts2);
       mrf.add_pairwise_factor(2,3,Potts2);
       mrf.add_pairwise_factor(0,3,Potts2);

       mrf.add_tightening_triplet(0,1,2);
       mrf.add_tightening_triplet(0,1,3);
       mrf.add_tightening_triplet(0,2,3);
       mrf.add_tightening_triplet(1,2,3);

       s.Solve();
       assert(false);
       test(s.lower_bound() == 1.0);
       }
     */
}
