#ifndef LPMP_MULTICUT_INSTANCE
#define LPMP_MULTICUT_INSTANCE

#include <array>
#include <vector>
#include <algorithm>
#include <cassert>
#include "union_find.hxx"
#include "config.hxx"

// TODO: split into .h and .cpp file

namespace LPMP {

   struct multicut_instance {
      enum class problem_type {multicut, correlation_clustering};

      struct weighted_edge : public std::array<std::size_t,2> { 
         weighted_edge()
            : std::array<std::size_t,2>({0,0}), cost(0.0) {}

         weighted_edge(const std::size_t i, const std::size_t j, const double c)
            : std::array<std::size_t,2>({i,j}), cost(c) {}

         double cost;
      };

      void add_edge(const std::size_t i, const std::size_t j, const double cost)
      {
         assert(i != j);
         no_nodes_ = std::max({no_nodes_,i+1,j+1});
         edges_.push_back(weighted_edge(std::min(i,j), std::max(i,j), cost));
      }

      const std::vector<weighted_edge>& edges() const { return edges_; }

      std::size_t no_nodes() const { return no_nodes_; }
      std::size_t no_edges() const { return edges_.size(); }

      void normalize()
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
                normalized_edges.back().cost += edges_[k].cost;
                ++k; 
             }
          }

          std::swap(normalized_edges, edges_); 
      }

      void transform_to_correlation_clustering()
      {
         if(problem_type_v == problem_type::correlation_clustering)
            return;
         problem_type_v = problem_type::correlation_clustering;

         involute_edges();
      }

      void transform_to_multicut()
      {
         if(problem_type_v == problem_type::multicut)
            return;
         problem_type_v = problem_type::multicut;

         involute_edges();
      }

      auto begin() { return edges_.begin(); }
      auto end() { return edges_.end(); }

      void add_to_constant(const double delta) { constant_ += delta; }

      class edge_labeling; // forward declaration

      class node_labeling : public std::vector<std::size_t> {
         public: 
         using std::vector<std::size_t>::vector;

         edge_labeling transform_to_edge_labeling(const multicut_instance& instance) const
         {
            edge_labeling output;
            output.transform_to_multicut();
            output.reserve(instance.no_edges());
            for(const auto& e : instance.edges()) {
               const bool cut = (*this)[e[0]] != (*this)[e[1]];
               output.push_back(cut); 
            }

            assert(std::abs(instance.evaluate(*this) == instance.evaluate(output)) <= eps);
            return output;
         } 
      };

      double evaluate(const node_labeling& l) const
      {
         assert(l.size() == no_nodes());
         double cost = constant_;
         if(problem_type_v == problem_type::multicut) {
            for(const auto& e : edges_) {
               if(l[e[0]] != l[e[1]])
                  cost += e.cost;
            }
         } else {
            for(const auto& e : edges_) {
               if(l[e[0]] == l[e[1]])
                  cost += e.cost;
            }
         }
         return cost;
      }

      class edge_labeling : public std::vector<char> {
         public:
            using std::vector<char>::vector;

            void transform_to_correlation_clustering()
            {
               if(problem_type_v == problem_type::correlation_clustering)
                  return;
               problem_type_v = problem_type::correlation_clustering;

               involute_edges();
            }

            bool check_primal_consistency(const multicut_instance& input) const
            {
               union_find uf(input.no_nodes());
               const char cut_value = problem_type_v == problem_type::multicut ? 1 : 0;
               for(std::size_t e=0; e<this->size(); ++e)
                  if((*this)[e] != cut_value)
                     uf.merge(input.edges()[e][0], input.edges()[e][1]);
               for(std::size_t e=0; e<this->size(); ++e)
                  if((*this)[e] == cut_value && uf.find(input.edges()[e][0]) == uf.find(input.edges()[e][1]))
                     return false;
               return true;
            }

            void transform_to_multicut()
            {
               if(problem_type_v == problem_type::multicut)
                  return;
               problem_type_v = problem_type::multicut;

               involute_edges();
            }

            node_labeling transform_to_node_labeling(const multicut_instance& instance) const
            {
               assert(this->size() == instance.no_edges());

               union_find uf(instance.no_nodes());
                  const char cut_value = problem_type_v == problem_type::multicut ? 1 : 0;

               for(std::size_t e=0; e<this->size(); ++e)
                     if((*this)[e] != cut_value)
                        uf.merge(instance.edges()[e][0], instance.edges()[e][1]); 

               node_labeling output(instance.no_nodes());
               for(std::size_t i=0; i<instance.no_nodes(); ++i)
                  output[i] = uf.find(i);

               assert(std::abs(instance.evaluate(*this) - instance.evaluate(output)) <= eps);
               return output;
            }

            problem_type get_problem_type() const { return problem_type_v; }

         private:
            problem_type problem_type_v = problem_type::multicut;

            void involute_edges()
            {
               for(auto& v : *this) {
                  assert(0 <= v && v <= 1);
                  v = 1 - v; 
               }
            }

      };

      double evaluate(const edge_labeling& l) const
      {
         assert(l.size() == no_edges());
         double cost = constant_;
         // TODO: check primal feasibility
         if(problem_type_v == l.get_problem_type()) {
            for(std::size_t e=0; e<no_edges(); ++e) {
               assert(l[e] == 0 || l[e] == 1);
               cost += this->edges()[e].cost * l[e];
            }
         } else {
            for(std::size_t e=0; e<no_edges(); ++e) {
               assert(l[e] == 0 || l[e] == 1);
               cost += this->edges()[e].cost * (1 - l[e]);
            } 
         }
         return cost; 
      }

      template<typename STREAM>
      void write_problem(STREAM& s) const
      {
         if(problem_type_v == problem_type::multicut)
            s << "MULTICUT\n";
         else if(problem_type_v == problem_type::correlation_clustering)
            s << "CORRELATION_CLUSTERING\n";

         for(const auto& e : edges_) {
            assert(e[0] < e[1]);
            s << e[0] << " " << e[1] << " " << e.cost << "\n"; 
         }
      }

      private:
      problem_type problem_type_v  = problem_type::multicut; 

      void involute_edges()
      {
         for(auto& e : edges_) {
            add_to_constant(e.cost);
            e.cost *= -1.0;
         } 
      }

      std::size_t no_nodes_ = 0;
      std::vector<weighted_edge> edges_; 

      double constant_ = 0.0;
   };

} // namespace LPMP

#endif // LPMP_MULTICUT_INSTANCE
