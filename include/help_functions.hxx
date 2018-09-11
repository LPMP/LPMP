#ifndef LPMP_HELP_FUNCTIONS_HXX
#define LPMP_HELP_FUNCTIONS_HXX

#include <string>
#include <map>
#include <set>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>
#include <assert.h>
#include <cstring>

#include <libgen.h>

namespace LPMP {

inline std::string LatexEscape(const std::string& s)
{
   std::string escaped = s;
   size_t start_pos = 0;
   while(start_pos < escaped.length() && (start_pos = escaped.find("_", start_pos)) != std::string::npos) {
      escaped.replace(start_pos, 1, "{\\_}");
      start_pos += 4;
   }
   return escaped;
}

// extract file name without extension, only works on linux for regular filenames.
inline std::string ExtractFilename(const std::string& path)
{
   char* pathTmp = new char[path.length()+1];
   std::strcpy(pathTmp,path.c_str());
   std::string file(basename(pathTmp));
   delete[] pathTmp;
   auto extensionPos = file.find_last_of("."); 
   if(extensionPos != std::string::npos) {
      file = file.substr(0,extensionPos);
   }
   return file;
}
//inline std::string ExtractFileName(const char* path)
//{
//   return ExtractFilename(std::string(path));
//}

inline int binary_compl(const int i)
{
   assert(i == 0 && i == 1);
   return 1-i;
}

template<class T>
int find_index(const T i, const std::vector<T>& vec)
{
   const int index = find(vec.begin(), vec.end(), i) - vec.begin();
   assert(index <= vec.size());
   return index;
}

template<class F, class I, typename ITERATOR>
void BuildIndexMaps(ITERATOR fIt, const ITERATOR fEndIt, std::map<F,I>& elemToIndex, std::map<I,F>& indexToElem)
{
   // do zrobienia: - reserve space, 
   //               - possibly use hash_map for speed
   for(std::size_t i=0; fIt+i!=fEndIt; ++i) {
      elemToIndex.insert(std::make_pair(*(fIt+i),i));
      indexToElem.insert(std::make_pair(i,*(fIt+i)));
   }
}

template<class T>
std::vector<T> GetSubVector(const std::vector<T>& v, const std::vector<size_t>& inds)
{
   std::vector<T> subVec(inds.size());
   for(size_t i=0; i<inds.size(); i++) {
      subVec[i] = v[inds[i]];
   }
   return subVec;
}

template<class VECTOR>
bool HasUniqueValues(const VECTOR& v)
{
   std::set<typename VECTOR::value_type> values;
   for(size_t i=0; i<v.size(); i++) {
      if(values.find( v[i] ) != values.end()) return false;
      else values.insert( v[i] );
   }
   return true;
}

template<typename T, typename ITERATOR>
T min_value(ITERATOR begin, ITERATOR end)
{
  T min_val = *begin;
  ++begin;
  for(; begin!=end; ++begin) {
    min_val = std::min(min_val, *begin);
  }
  return min_val;
}

// find out smallest and second smallest without branching
template<typename T, typename ITERATOR>
std::array<T,2> two_smallest_elements(ITERATOR begin, ITERATOR end)
{
  T smallest = std::numeric_limits<T>::infinity();
  T second_smallest = std::numeric_limits<T>::infinity();
  for(; begin!=end; ++begin) {
    const auto min = std::min(smallest, *begin);
    const auto max = std::max(smallest, *begin);
    smallest = min;
    second_smallest = std::min(max, second_smallest);
  }
  return {smallest, second_smallest};
}

template<typename T, typename ITERATOR>
std::array<T,2> two_largest_elements(ITERATOR begin, ITERATOR end)
{
  T largest = -std::numeric_limits<T>::infinity();
  T second_largest = -std::numeric_limits<T>::infinity();
  for(; begin!=end; ++begin) {
    const auto max = std::max(largest, *begin);
    const auto min = std::min(largest, *begin);
    largest = max;
    second_largest = std::max(min, second_largest);
  }
  return {largest, second_largest};
}

// return indices belonging to the three smallest entries
template<class T>
std::tuple<size_t,size_t,size_t> MinThreeIndices(const std::vector<T>& v)
{
   assert(v.size() > 2);
   std::vector<size_t> ind(v.size());
   for(size_t i=0; i<ind.size(); i++)  ind[i] = i;
   // do zrobienia: use partial_sort here
   std::sort(ind.begin(), ind.end(), [&](const size_t i, const size_t j)->bool {return v[i] < v[j];} );
   return std::tuple<size_t,size_t,size_t>(ind[0], ind[1], ind[2]);
}

template<class T, class A>
std::pair<T,T> SmallestValues(const A& v)
{
   assert(v.size() > 1);
   T min_val = std::numeric_limits<T>::infinity();
   T second_min_val = std::numeric_limits<T>::infinity();
   for(std::size_t i=0; i<v.size(); ++i) {
      
      assert(false); 
   }
   return std::make_pair(min_val, second_min_val);
}


template<class T>
void NormalizeVector(std::vector<T>& v)
{
   //const T mean = std::accumulate(v.begin(), v.end(),0)/v.size();
   return; // do zrobienia: why is it not applicable in marg_message?
   const T mean = *std::min_element(v.begin(), v.end());
   for(size_t i=0; i<v.size(); i++) {
      v[i] -= mean;
   }
}

// compute permutation of v1 onto v2
template<typename VECTOR>
std::vector<std::size_t> ComputePermutation(const VECTOR& v1, const VECTOR& v2)
{
   assert(v1.size() == v2.size());
   //assert(std::is_unique(v1));
   std::map<typename std::remove_reference<typename std::remove_cv<decltype(v1[0])>::type>::type,std::size_t> m;
   for(std::size_t i=0; i<v1.size(); ++i) {
      m.insert(std::make_pair(v1[i],i));
   }
   std::vector<std::size_t> perm(v1.size());
   for(std::size_t i=0; i<v2.size(); ++i) {
      perm[i] = m[v2[i]]; 
   }
   return perm;
}


template<
class InputIt1, class InputIt2, class OutputIt, class Compare, class Merge >
   static OutputIt set_intersection_merge
(
 InputIt1 first1, InputIt1 last1,
 InputIt2 first2, InputIt2 last2,
 OutputIt d_first, Compare comp, Merge merge
 )
{
   while (first1 != last1 && first2 != last2)
   {
      if (comp(*first1, *first2))
         ++first1;
      else
      {
         if (!comp(*first2, *first1))
            *d_first++ = merge(*first1++, *first2);
         ++first2;
      }
   }
   return d_first;
}




} // end namespace LPMP

#endif // LPMP_HELP_FUNCTIONS_HXX
