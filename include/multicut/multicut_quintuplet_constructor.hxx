#pragma once

#include "cut_base/cut_base_quintuplet_constructor.hxx"
#include "multicut_odd_bicycle_wheel_packing.h"

namespace LPMP {

template<typename QUADRUPLET_CONSTRUCTOR, 
   typename QUINTUPLET_FACTOR,
   typename QUADRUPLET_QUINTUPLET_MESSAGE_0123,
   typename QUADRUPLET_QUINTUPLET_MESSAGE_0124,
   typename QUADRUPLET_QUINTUPLET_MESSAGE_0134,
   typename QUADRUPLET_QUINTUPLET_MESSAGE_0234,
   typename QUADRUPLET_QUINTUPLET_MESSAGE_1234
   >
class multicut_quintuplet_constructor : public cut_base_quadruplet_quintuplet_constructor<
					       //multicut_quintuplet_constructor< QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>,
					       QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>
{
public:
   using FMC = typename QUADRUPLET_CONSTRUCTOR::FMC;
   using base_constructor = cut_base_quadruplet_quintuplet_constructor<QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>;

   using quintuplet_factor_container = QUINTUPLET_FACTOR;
   using odd_3_wheel_quintuplet_message_0123_container = QUADRUPLET_QUINTUPLET_MESSAGE_0123;
   using odd_3_wheel_quintuplet_message_0124_container = QUADRUPLET_QUINTUPLET_MESSAGE_0124;
   using odd_3_wheel_quintuplet_message_0134_container = QUADRUPLET_QUINTUPLET_MESSAGE_0134;
   using odd_3_wheel_quintuplet_message_0234_container = QUADRUPLET_QUINTUPLET_MESSAGE_0234;
   using odd_3_wheel_quintuplet_message_1234_container = QUADRUPLET_QUINTUPLET_MESSAGE_1234;

   template<typename SOLVER>
   multicut_quintuplet_constructor(SOLVER& s) : base_constructor(s) {}

   std::size_t Tighten(const std::size_t no_constraints_to_add);
   void ComputePrimal();

   template<typename ITERATOR>
   auto* add_higher_order_quintuplet(
         const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3, const std::size_t i4,
         ITERATOR cost_begin, ITERATOR cost_end)
   {
      auto* f = this->add_quintuplet(i0,i1,i2,i3);
      auto* t = f->get_factor();
      assert(std::distance(cost_begin, cost_end) + 1 == t->size());
      this->add_to_constant(*cost_begin);

      std::size_t i = 0;
      for(auto it = cost_begin+1; it != cost_end; ++it, ++i) {
         (*t)[i] = *it - *cost_begin; 
      }

      return f;
   }

   // ij is the axle, uv is the wheel edge
   // labelings with edge ij and uv cut and (i) iu and jv or (ii) iv and ju not cut must have smaller cost than other configurations
   double compute_edge_cost(const std::size_t i, const std::size_t j, const std::size_t u, const std::size_t v)
   {
      std::array<std::size_t,4> idx{i,j,u,v};
      std::sort(idx.begin(), idx.end());
      // can it happen that odd bicycle wheel inequalities can tighten the polytope but odd wheel inequalities are not also violated? Only in this case the below construction makes sense
      if(!this->has_odd_3_wheel_factor(idx[0], idx[1], idx[2], idx[3])) { // create a fake odd 3 wheel factor and reparamaetrize all underlying triplets into it

	      multicut_odd_3_wheel_factor f;

	      // possibly do not create new messages, but use static functions in these messages directly
	      if(this->HasTripletFactor(idx[0],idx[1],idx[2])) {
		      const auto& t = *(this->GetTripletFactor(idx[0], idx[1], idx[2])->get_factor());
		      multicut_triplet_odd_3_wheel_message_012 m;
		      m.RepamRight(f,t);
	      }

	      if(this->HasTripletFactor(idx[0],idx[1],idx[3])) {
		      const auto& t = *(this->GetTripletFactor(idx[0], idx[1], idx[3])->get_factor());
		      multicut_triplet_odd_3_wheel_message_013 m;
		      m.RepamRight(f,t);
	      }

	      if(this->HasTripletFactor(idx[0],idx[2],idx[3])) {
		      const auto& t = *(this->GetTripletFactor(idx[0], idx[2], idx[3])->get_factor());
		      multicut_triplet_odd_3_wheel_message_023 m;
		      m.RepamRight(f,t);
	      }

            if(this->HasTripletFactor(idx[1],idx[2],idx[3])) {
               const auto& t = *(this->GetTripletFactor(idx[1], idx[2], idx[3])->get_factor());
               multicut_triplet_odd_3_wheel_message_123 m;
               m.RepamRight(f,t);
            }
            return this->compute_edge_cost_from_odd_3_wheel(i,j,u,v, f);

      } else {
         auto& f = *(this->get_odd_3_wheel_factor(idx[0], idx[1], idx[2], idx[3])->get_factor());
         return this->compute_edge_cost_from_odd_3_wheel(i,j,u,v, f);
      }
   }

   double compute_edge_cost_from_odd_3_wheel(const std::size_t i, const std::size_t j, const std::size_t u, const std::size_t v, const multicut_odd_3_wheel_factor& f)
   {
      assert(i < j && u < v);
      std::array<std::size_t,4> idx{i,j,u,v};
      std::sort(idx.begin(), idx.end());
      double min_participating_labelings = std::numeric_limits<double>::infinity();
      double min_non_participating_labelings;
      // non participating labelings have always != 4 cut edges
      min_non_participating_labelings = std::min({0.0, f[0], f[1], f[2], f[3], f[7], f[8], f[9], f[10], f[11], f[12], f[13]});
      // for participating edges, the axle and wheel edge must be cut
      if(idx[0] == i && idx[1] == j) { // first and sixth edge
         assert(idx[2] == u && idx[3] == v);
         min_participating_labelings = std::min(f[5], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[4]);
      } else if(idx[0] == i && idx[2] == j) { // second and fifth
         assert(idx[1] == u && idx[3] == v);
         min_participating_labelings = std::min(f[4], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[5]);
      } else if(idx[0] == i && idx[3] == j) { // third and fourth
         assert(idx[1] == u && idx[2] == v);
         min_participating_labelings = std::min(f[4], f[5]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[6]);
      } else if(idx[1] == i && idx[2] == j) { // fourth and third
         assert(idx[0] == u && idx[3] == v);
         min_participating_labelings = std::min(f[4], f[5]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[6]);
      } else if(idx[1] == i && idx[3] == j) { // fifth and second
         assert(idx[0] == u && idx[2] == v);
         min_participating_labelings = std::min(f[4], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[5]);
      } else if(idx[2] == i && idx[3] == j) { // sixth and first
         assert(idx[0] == u && idx[1] == v);
         min_participating_labelings = std::min(f[5], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[4]);
      } else {
         assert(false);
      }
      //assert(false);

      return min_participating_labelings - min_non_participating_labelings; 
   }

   /*
      static std::array<double,2> odd_3_wheel_xxx_to_do(const std::size_t i, const std::size_t j, const std::size_t u, const std::size_t v, const typename base_constructor::odd_3_wheel_factor& f)
      {
      double min_participating_labelings = std::numeric_limits<double>::infinity();
      double min_non_participating_labelings;

      min_non_participating_labelings = std::min({0.0, f[0], f[1], f[2], f[3], f[7], f[8], f[9], f[10], f[11], f[12], f[13]});

      // for participating edges, the axle and wheel edge must be cut
      if(idx[0] == i && idx[1] == j) { // first and sixth edge
         assert(idx[2] == u && idx[3] == v);
         min_participating_labelings = std::min(f[5], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[4]);
      } else if(idx[0] == i && idx[2] == j) { // second and fifth
         assert(idx[1] == u && idx[3] == v);
         min_participating_labelings = std::min(f[4], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[5]);
      } else if(idx[0] == i && idx[3] == j) { // third and fourth
         assert(idx[1] == u && idx[2] == v);
         min_participating_labelings = std::min(f[4], f[5]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[6]);
      } else if(idx[1] == i && idx[2] == j) { // fourth and third
         assert(idx[0] == u && idx[3] == v);
         min_participating_labelings = std::min(f[4], f[5]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[6]);
      } else if(idx[1] == i && idx[3] == j) { // fifth and second
         assert(idx[0] == u && idx[2] == v);
         min_participating_labelings = std::min(f[4], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[5]);
      } else if(idx[2] == i && idx[3] == j) { // sixth and first
         assert(idx[0] == u && idx[1] == v);
         min_participating_labelings = std::min(f[5], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[4]);
      } else {
         assert(false);
      }

      return {min_participating_labelings, min_non_participating_labelings};
   }
*/

   using triangle_intersection_type = std::tuple<std::size_t,std::size_t, typename base_constructor::triplet_factor*, typename base_constructor::triplet_factor*>;

   void compute_triangles( // triangle are pyramids, though, search for better name
         const std::size_t i, const std::size_t j, const double minTh, 
         const typename base_constructor::triplet_connections& connected_triplets,
         std::vector<triangle_intersection_type>& common_edges,
         std::unordered_map<std::size_t,std::size_t>& origToCompressedNode, 
         std::vector<std::size_t>& compressedToOrigNode, 
         std::vector<std::tuple<std::size_t,std::size_t,double>>& compressedEdges)
   {
      assert(i < j);
      origToCompressedNode.clear();
      compressedToOrigNode.clear();
      compressedEdges.clear();

      auto merge = [](const auto a, const auto b) -> triangle_intersection_type { 
         assert(a.nodes[0] == b.nodes[0] && a.nodes[1] == b.nodes[1]);
         return std::make_tuple(a.nodes[0], a.nodes[1], a.f, b.f); 
      };

      // find all triangles ijk

      // find all edges uv such that there exist edge triplets iuv and juv. 
      // this is done by sorting all triplets which have node i and node j, and intersecting the set
      auto intersects_iter_end = set_intersection_merge(
            connected_triplets[i].begin(), connected_triplets[i].end(),
            connected_triplets[j].begin(), connected_triplets[j].end(),
            common_edges.begin(), [](const auto& a, const auto& b) { return a.operator<(b); }, merge);

      for(auto n=common_edges.begin(); n != intersects_iter_end; ++n) {
         const std::size_t u = std::get<0>(*n);
         const std::size_t v = std::get<1>(*n);
         assert(u < v);
         const auto& iuv = std::get<2>(*n)->get_factor();
         const auto& juv = std::get<3>(*n)->get_factor();

         const double dual_increase = compute_edge_cost(i,j,u,v);

         if(dual_increase >= minTh) { // add edge uv to bipartite graph

            if(origToCompressedNode.find(u) == origToCompressedNode.end()) {
               origToCompressedNode.insert(std::make_pair(u, origToCompressedNode.size()));
               compressedToOrigNode.push_back(u);
            }
            if(origToCompressedNode.find(v) == origToCompressedNode.end()) {
               origToCompressedNode.insert(std::make_pair(v, origToCompressedNode.size()));
               compressedToOrigNode.push_back(v);
            }
            const std::size_t uc = origToCompressedNode[u];
            const std::size_t vc = origToCompressedNode[v];
            assert(uc != vc);
            compressedEdges.push_back(std::make_tuple(uc,vc, dual_increase));
         }
      } 
   }

   /*
   template<typename ITERATOR>
   void triangulate_quintuplet(const std::size_t i, const std::size_t j, const double cost, ITERATOR path_begin, ITERATOR path_end, std::vector<bicycle_candidate>& candidates)
   {
      assert(i < j);
      assert(std::distance(path_begin, path_end) >= 3);
      this->cycle_normal_form(path_begin, path_end);
      const std::size_t first_node = *path_begin;
      for(auto it=path_begin+1; it+1!=path_end; ++it) {
         std::array<std::size_t,5> nodes({i,j,first_node, *it, *(it+1)});
         std::sort(nodes.begin(), nodes.end());
         assert(HasUniqueValues(nodes));
         candidates.push_back({nodes, cost});
      }
   }
   */


   /*
   std::size_t find_violated_quintuplets(const std::size_t max_factors_to_add)
   {
      if(this->number_of_edges() > 2) {
         // preprocessing: sort triplets for fast intersection later
         auto connected_triplets = this->compute_connected_triplets();

         std::vector<std::tuple<std::size_t,double>> threshold(this->unaryFactorsVector_.size()); // edge number and threshold
         std::vector<bicycle_candidate> odd_bicycle_candidates;
         // given a cut axle edge (negative cost), find cut wheel edges such that among the four spokes exactly two are cut and two are connected.
     
//#pragma omp parallel
         {
            using intersection_type = std::tuple<std::size_t,std::size_t, typename base_constructor::triplet_factor*, typename base_constructor::triplet_factor*>;
            std::vector<intersection_type> common_edges(this->number_of_edges()); // possibly this is a bit large!
            auto merge = [](const auto a, const auto b) -> intersection_type { 
               assert(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) == std::get<1>(b));
               return std::make_tuple(std::get<0>(a), std::get<1>(a), std::get<2>(a), std::get<2>(b)); 
            };
            std::vector<bicycle_candidate> odd_bicycle_candidates_local;

            std::unordered_map<std::size_t,std::size_t> origToCompressedNode;
            std::vector<std::size_t> compressedToOrigNode;
            std::vector<std::tuple<std::size_t,std::size_t,double>> compressedEdges;

//#pragma omp for schedule(guided) nowait
            for(std::size_t e=0; e<this->unaryFactorsVector_.size(); ++e) {

               // edge ij will be treated as axle of odd bicycle wheel
               const std::size_t i = std::get<0>(this->unaryFactorsVector_[e])[0];
               const std::size_t j = std::get<0>(this->unaryFactorsVector_[e])[1];
               const double cost_ij = std::get<1>(this->unaryFactorsVector_[e])->get_factor()->operator[](0); 
               if(cost_ij < -eps) {
                  //origToCompressedNode.clear();
                  //compressedToOrigNode.clear();
                  //compressedEdges.clear(); 

                  compute_triangles(i, j, eps, connected_triplets, common_edges, origToCompressedNode, compressedToOrigNode, compressedEdges); 
                  const double th = this->compute_odd_cycle_threshold(origToCompressedNode, compressedToOrigNode, compressedEdges);

                  if(th > eps) {
                     auto path = this->compute_path_in_bipartite_graph(origToCompressedNode, compressedToOrigNode, compressedEdges, th);
                     //triangulate_quintuplet(i,j, th, path.begin(), path.end(), odd_bicycle_candidates_local);
                  }
               } 
            }
//#pragma omp critical
            odd_bicycle_candidates.insert(odd_bicycle_candidates.end(), odd_bicycle_candidates_local.begin(), odd_bicycle_candidates_local.end());
         }

         std::sort(odd_bicycle_candidates.begin(), odd_bicycle_candidates.end(), [](const auto& a, const auto& b) { return a.cost > b.cost; });
         std::size_t no_factors_added = 0;
         for(std::size_t i=0; i<odd_bicycle_candidates.size(); ++i) {
            if(!this->has_quintuplet( odd_bicycle_candidates[i].idx )) {
               this->add_quintuplet( odd_bicycle_candidates[i].idx );
               ++no_factors_added;
               if(no_factors_added >= max_factors_to_add) {
                  break;
               }
            } 
         }
         return no_factors_added;
      } else {
         return 0;
      }
   }
   */

private:

   std::unordered_map<std::array<std::size_t,5>, quintuplet_factor_container*> quintuplet_factors_;

   //std::unordered_map<std::vector<std::tuple<std::size_t,std::size_t, typename base_constructor::odd_3_wheel_factor_container*>>> odd_3_wheel_factor_by_indices_; // if odd 3 wheel factor with indices (i1,i2,i3,i4) exists, then (i1,i2,i3,i4) will be in the hash indexed by all two-subsets of the indices.

};

    template<typename QUADRUPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename QUADRUPLET_QUINTUPLET_MESSAGE_0123, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234>
std::size_t multicut_quintuplet_constructor<QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>::Tighten(const std::size_t no_constraints_to_add)
{
    this->send_messages_to_quadruplets();
    quadruplet_multicut_instance mc = this->template export_quadruplets<quadruplet_multicut_instance>();
    const odd_bicycle_wheel_packing obwp = compute_multicut_odd_bicycle_wheel_packing(mc);
    if(debug())
        std::cout << "found " << obwp.no_odd_bicycle_wheels() << " odd bicycle wheels\n";
    for(std::size_t c=0; c<obwp.no_odd_bicycle_wheels(); ++c) {
        const auto [cycle_begin, cycle_end] = obwp.get_cycle(c);
        const double odd_bicycle_wheel_weight = obwp.get_odd_bicycle_wheel_weight(c);
        const auto axle = obwp.get_axle(c);
        this->triangulate_odd_bicycle_wheel(axle, cycle_begin, cycle_end, odd_bicycle_wheel_weight);
    }

    const std::size_t base_constraints_added = base_constructor::Tighten(no_constraints_to_add);

    return base_constraints_added + obwp.no_odd_bicycle_wheels();
}
    template<typename QUADRUPLET_CONSTRUCTOR, typename QUINTUPLET_FACTOR, typename QUADRUPLET_QUINTUPLET_MESSAGE_0123, typename QUADRUPLET_QUINTUPLET_MESSAGE_0124, typename QUADRUPLET_QUINTUPLET_MESSAGE_0134, typename QUADRUPLET_QUINTUPLET_MESSAGE_0234, typename QUADRUPLET_QUINTUPLET_MESSAGE_1234>
void multicut_quintuplet_constructor<QUADRUPLET_CONSTRUCTOR, QUINTUPLET_FACTOR, QUADRUPLET_QUINTUPLET_MESSAGE_0123, QUADRUPLET_QUINTUPLET_MESSAGE_0124, QUADRUPLET_QUINTUPLET_MESSAGE_0134, QUADRUPLET_QUINTUPLET_MESSAGE_0234, QUADRUPLET_QUINTUPLET_MESSAGE_1234>::ComputePrimal()
{
    this->send_messages_to_quadruplets();
    base_constructor::ComputePrimal(); 
}

} // namespace LPMP
