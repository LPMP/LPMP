#include "test.h"
#include "topological_sort.hxx"
#include <random>
#include <algorithm>

using namespace LPMP;

void test_topological_sort(const std::size_t n, const std::size_t m, std::random_device& rd)
{
   // generate random permutation of [1,n]
   std::vector<std::size_t> permutation(n);
   std::iota(permutation.begin(), permutation.end(), 0);
   std::mt19937 gen(rd());
   std::shuffle(permutation.begin(), permutation.end(), gen);

   // sample randomly m pairs of indices in permutation and add to topological sorting graph
   std::bernoulli_distribution dist{double(m)/double(0.5*(n*(n-1)))};
   Topological_Sort::Graph tsg(n);
   std::vector<std::array<std::size_t,2>> precedence_relations; 
   for(std::size_t i=0; i<n; ++i) {
      for(std::size_t j=i+1; j<n; ++j) {
         if(dist(gen)) {
            precedence_relations.push_back({permutation[i],permutation[j]});
            tsg.addEdge(permutation[i], permutation[j]);
         }
      }
   }
   
   auto sorting = tsg.topologicalSort();

   // check whether all precedence relations are fulfilled
   test(tsg.sorting_valid(sorting));

   tsg.addEdge(0,n/2);
   tsg.addEdge(n/2,n-1);
   tsg.addEdge(n-1,0);
   try { // must throw for the test to succeed
      tsg.topologicalSort();
      throw test_exception("graph not DAG not recognized."); // non-standard exception, will not be caught
   } catch (std::exception& error) {} 
}

int main(int argc, char**argv)
{
    std::random_device rd{};
    for(std::size_t n=5; n<100; n+= 5) {
       for(std::size_t m=n/4; m<n*(n-1)/4; m+=n/5) {
          test_topological_sort(n,m,rd);
       }
    } 
}
