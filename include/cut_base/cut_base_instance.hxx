#pragma once 

#include <array>
#include <vector>
#include <cassert>
#include <algorithm>
#include <stack>
#include "vector.hxx"
#include "graph.hxx"
#include "union_find.hxx"

namespace LPMP {

    class cut_base_node_labeling; // forward declaration
    class cut_base_edge_labeling; // forward declaration

    struct cut_base_instance {

        struct edge_factor : public array<double,1> {
            operator double() const { return (*this)[0]; }
        };

        struct weighted_edge : public std::array<std::size_t,2> { 
            weighted_edge()
                : std::array<std::size_t,2>({0,0}), cost({0.0}) {}

            weighted_edge(const std::size_t i, const std::size_t j, const double c)
                : std::array<std::size_t,2>({i,j}), cost({c}) {}

            edge_factor cost;
        };

        void add_edge(const std::size_t i, const std::size_t j, const double cost);
        const std::vector<weighted_edge>& edges() const { return edges_; }
        std::vector<weighted_edge>& edges() { return edges_; }
        std::size_t no_nodes() const { return no_nodes_; }
        std::size_t no_edges() const { return edges_.size(); }

        void normalize(); // merge parallel edges
        bool normalized() const; // check if there are parallel edges

        auto begin() { return edges_.begin(); }
        auto end() { return edges_.end(); }

        void add_to_constant(const double delta) { constant_ += delta; }
        double constant() const { return constant_; }

        double evaluate(const cut_base_node_labeling& l) const; 
        double evaluate(const cut_base_edge_labeling& l) const;
        double lower_bound() const;

        bool graph_connected() const;
        void shift_to_zero_offset();

        template<typename STREAM>
            void write_problem(STREAM& s) const;

        protected:
        std::size_t no_nodes_ = 0;
        std::size_t min_node_ = std::numeric_limits<std::size_t>::max();
        std::vector<weighted_edge> edges_; 

        double constant_ = 0.0;
    };

    class cut_base_node_labeling : public std::vector<std::size_t> {
        public: 
            using std::vector<std::size_t>::vector;
    };

    class cut_base_edge_labeling : public std::vector<char> {
        public:
            using std::vector<char>::vector;

            cut_base_edge_labeling(const cut_base_instance& instance, const cut_base_node_labeling& labeling);
            cut_base_edge_labeling(const cut_base_instance& instance, union_find& uf);
    };


    template<typename TRIPLET_FACTOR>
        class triplet_cut_base_instance : public cut_base_instance {
            public:
                using triplet_factor = TRIPLET_FACTOR;
                // TODO: or store below the edges which the triplet covers?
                struct cut_base_triplet : public std::array<std::size_t,3> { // the three nodes of the triplet
                    triplet_factor cost;
                };

                void add_triplet(std::array<std::size_t,3> nodes, triplet_factor t);
                double lower_bound() const; 
                const std::vector<cut_base_triplet>& triplets() const { return triplets_; }
                std::vector<cut_base_triplet>& triplets() { return triplets_; }

                bool normalized() const;

            private:
                std::vector<cut_base_triplet> triplets_;
        };

    template<typename QUADRUPLET_FACTOR, typename BASE_INSTANCE>
        class quadruplet_cut_base_instance : public BASE_INSTANCE {
            public:
                using quadruplet_factor = QUADRUPLET_FACTOR;
                // TODO: or store below the edges which the triplet covers?
                struct cut_base_quadruplet : public std::array<std::size_t,4> { // the three nodes of the quadruplet
                    quadruplet_factor cost;
                };

                void add_quadruplet(std::array<std::size_t,4> nodes, const QUADRUPLET_FACTOR& f);
                double lower_bound() const; 
                const std::vector<cut_base_quadruplet>& quadruplets() const { return quadruplets_; }
                std::vector<cut_base_quadruplet>& quadruplets() { return quadruplets_; }

                bool normalized() const;
            private: 
                std::vector<cut_base_quadruplet> quadruplets_; 
        };


    // implementation

    inline void cut_base_instance::add_edge(const std::size_t i, const std::size_t j, const double cost)
    {
        assert(i != j);
        no_nodes_ = std::max({no_nodes_,i+1,j+1});
        min_node_ = std::min({i, j, min_node_});
        edges_.push_back(weighted_edge(std::min(i,j), std::max(i,j), cost));
    }

    inline void cut_base_instance::normalize()
    {
        // merge parallel edges
        std::sort(edges_.begin(), edges_.end(), [](const auto& e1, const auto& e2) { 
                return std::lexicographical_compare(e1.begin(), e1.end(), e2.begin(), e2.end());
                });

        std::vector<weighted_edge> normalized_edges; 

        // merge matching edge copies and add them to edges
        for(std::size_t k=0; k<edges_.size();) {
            normalized_edges.push_back(edges_[k]);
            ++k; 
            while(k<edges_.size() && normalized_edges.back()[0] == edges_[k][0] && normalized_edges.back()[1] == edges_[k][1]) {
                normalized_edges.back().cost[0] += edges_[k].cost[0];
                ++k; 
            }
        }

        std::swap(normalized_edges, edges_); 
        assert(normalized());
    }

    inline bool cut_base_instance::normalized() const
    {
        auto edges_copy = edges_;
        std::sort(edges_copy.begin(), edges_copy.end(), [](const auto& e1, const auto& e2) { 
            return std::lexicographical_compare(e1.begin(), e1.end(), e2.begin(), e2.end());
        });
        return std::unique(edges_copy.begin(), edges_copy.end(), [](const auto& e1, const auto& e2) {  
                return std::equal(e1.begin(), e1.end(), e2.begin());
                }) == edges_copy.end(); 
    }

    template<typename STREAM>
        void cut_base_instance::write_problem(STREAM& s) const
        {
            for(const auto& e : edges_)
                s << e[0] << " " << e[1] << " " << e.cost << "\n";
        }

    inline cut_base_edge_labeling::cut_base_edge_labeling(const cut_base_instance& instance, union_find& uf)
    {
        assert(instance.no_nodes() == uf.size());
        this->reserve(instance.no_edges());
        for(const auto& e : instance.edges()) {
            if(uf.connected(e[0], e[1]))
                this->push_back(0);
            else
                this->push_back(1);
        }
    }

    inline cut_base_edge_labeling::cut_base_edge_labeling(const cut_base_instance& instance, const cut_base_node_labeling& labeling)
    {
        this->reserve(instance.no_edges());
        for(const auto& e : instance.edges()) {
            const bool cut = labeling[e[0]] != labeling[e[1]];
            this->push_back(cut); 
        }

        assert(std::abs(instance.evaluate(*this) - instance.evaluate(labeling)) <= 1e-8);
    } 

    inline double cut_base_instance::evaluate(const cut_base_node_labeling& l) const
    {
        assert(l.size() == no_nodes());
        double cost = constant_;
        for(const auto& e : edges_) {
            if(l[e[0]] != l[e[1]])
                cost += e.cost;
        }
        return cost;
    }

    inline double cut_base_instance::lower_bound() const
    {
        double lb = constant_;
        for(const auto& e : edges_) {
            lb += std::min(e.cost[0], 0.0);
        }
        return lb;
    }

    inline bool cut_base_instance::graph_connected() const
    {
        union_find uf(no_nodes());
        for(const auto e : edges())
            uf.merge(e[0], e[1]);
        return uf.count() == 1;
    }

    inline void cut_base_instance::shift_to_zero_offset()
    {
// TODO: this needs to be implemented for triplets, quadruplets and quintuplets as well!

        if(min_node_ > 0) {
            for(auto& e : edges_) {
                e[0] -= min_node_;
                e[1] -= min_node_;
            }
        }
        no_nodes_ -= min_node_;
        min_node_ = 0;
    }


    inline double cut_base_instance::evaluate(const cut_base_edge_labeling& l) const
    {
        assert(l.size() == no_edges());
        double cost = constant_;
        // TODO: check primal feasibility
        for(std::size_t e=0; e<no_edges(); ++e) {
            assert(l[e] == 0 || l[e] == 1);
            cost += this->edges()[e].cost * l[e];
        }
        return cost; 
    }

    /* ---------------
       lifted multicut
       ---------------- */

    /*
       bool lifted_cut_base_instance::edge_labeling::check_primal_consistency(const lifted_cut_base_instance& input) const
       {
    //collect connectivity information with union find w.r.t. base edges
    union_find uf(input.base.no_nodes());

    for(std::size_t e=0; e<base.size(); ++e)
    if(base[e] == 0)
    uf.merge(input.base.edges()[e][0],input.base.edges()[e][1]);

    for(std::size_t e=0; e<lifted.size(); ++e)
    if(lifted[e] == 0 && uf.connected(input.lifted.edges()[e][0], input.lifted.edges()[e][1]))
    return false;

    // check feasibility of multicut w.r.t. base and lifted edges
    for(std::size_t e=0; e<lifted.size(); ++e)
    if(lifted[e] == 0)
    uf.merge(input.lifted.edges()[e][0],input.lifted.edges()[e][1]);

    for(std::size_t e=0; e<base.size(); ++e)
    if(base[e] == 1 && uf.find(input.base.edges()[e][0]) == uf.find(input.base.edges()[e][1]))
    return false;

    for(std::size_t e=0; e<base.size(); ++e)
    if(lifted[e] == 1 && uf.find(input.lifted.edges()[e][0]) == uf.find(input.lifted.edges()[e][1]))
    return false;

    return true;
    }
    */

    template<typename TRIPLET_FACTOR>
        double triplet_cut_base_instance<TRIPLET_FACTOR>::lower_bound() const
        {
            double lb = cut_base_instance::lower_bound();
            for(const auto& t : triplets_) {
                lb += t.cost.LowerBound();
            }
            return lb;
        }

    template<typename TRIPLET_FACTOR>
        void triplet_cut_base_instance<TRIPLET_FACTOR>::add_triplet(std::array<std::size_t,3> nodes, TRIPLET_FACTOR t)
        {
            assert(nodes[0] < nodes[1] && nodes[1] < nodes[2] && nodes[2] <= this->no_nodes());
            triplets_.push_back({nodes, t});
        }

    template<typename TRIPLET_FACTOR>
        bool triplet_cut_base_instance<TRIPLET_FACTOR>::normalized() const
        {
            if(!cut_base_instance::normalized()) 
                return false;

            auto triplets_copy = triplets_;
            std::sort(triplets_copy.begin(), triplets_copy.end(), [](const auto& t1, const auto& t2) {
                    return std::lexicographical_compare(t1.begin(), t1.end(), t2.begin(), t2.end());
                    });

            return std::unique(triplets_copy.begin(), triplets_copy.end(), [](const auto& t1, const auto& t2) {
                    return std::equal(t1.begin(), t1.end(), t2.begin());
                    }) == triplets_copy.end();
        }

    template<typename QUADRUPLET_FACTOR, typename BASE_INSTANCE>
        double quadruplet_cut_base_instance<QUADRUPLET_FACTOR, BASE_INSTANCE>::lower_bound() const
        {
            double lb = BASE_INSTANCE::lower_bound();
            for(const auto& t : quadruplets_) {
                lb += t.cost.LowerBound();
            }
            return lb;
        }

    template<typename QUADRUPLET_FACTOR, typename BASE_INSTANCE>
        void quadruplet_cut_base_instance<QUADRUPLET_FACTOR, BASE_INSTANCE>::add_quadruplet(std::array<std::size_t,4> nodes, const QUADRUPLET_FACTOR& f)
        {
            assert(nodes[0] < nodes[1] && nodes[1] < nodes[2] && nodes[2] < nodes[3] && nodes[3] <= this->no_nodes());
            quadruplets_.push_back({nodes, f});
        }

    template<typename QUADRUPLET_FACTOR, typename BASE_INSTANCE>
        bool quadruplet_cut_base_instance<QUADRUPLET_FACTOR, BASE_INSTANCE>::normalized() const
        {
            if(!BASE_INSTANCE::normalized()) 
                return false;

            auto quadruplets_copy = quadruplets_;
            std::sort(quadruplets_copy.begin(), quadruplets_copy.end(), [](const auto& q1, const auto& q2) {
                    return std::lexicographical_compare(q1.begin(), q1.end(), q2.begin(), q2.end());
                    });

            return std::unique(quadruplets_copy.begin(), quadruplets_copy.end(), [](const auto& q1, const auto& q2) {
                    return std::equal(q1.begin(), q1.end(), q2.begin());
                    }) == quadruplets_copy.end();
        }

} // namespace LPMP
