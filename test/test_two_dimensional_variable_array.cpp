#include "test.h"
#include "two_dimensional_variable_array.hxx"
#include <random>

using namespace LPMP;

template<typename T>
void test_two_dimensional_variable_array(T val_1, T val_2)
{
   // random vector tests
   std::random_device rd{};

   for(std::size_t n=2; n<100; ++n) {
      std::mt19937 gen{rd()};
      std::uniform_int_distribution<> dist{0,1000};
      std::vector<std::size_t> size(n);
      for(auto& x : size) {
         x = dist(gen);
      }

      two_dim_variable_array<T> array(size.begin(), size.begin() + size.size()/2, val_1);
      test(array.size() == size.size()/2);
      {
         auto it = array.begin();
         for(std::size_t i=0; i<array.size(); ++i, ++it) {
            test(array[i].size() == size[i]); 
            test((*it).size() == size[i]); 
         }
         test(it == array.end());
      }
      {
         auto it = array.begin();
         for(std::size_t i=0; i<array.size(); ++i, ++it) { 
            for(std::size_t j=0; j<array[i].size(); ++j) { 
               test(array(i,j) == val_1);
               test(array[i][j] == val_1);
               test((*it)[j] == val_1);
            }
         }
      }

      array.resize(size.begin(), size.end(), 2.0);
      test(array.size() == size.size());
      for(std::size_t i=0; i<array.size(); ++i) { test(array[i].size() == size[i]); }
      for(std::size_t i=0; i<array.size(); ++i) { 
         for(std::size_t j=0; j<array[i].size(); ++j) { 
            if(i < size.size()/2) {
               test(array(i,j) == val_1);
               test(array[i][j] == val_1);
            } else {
               test(array(i,j) == val_2);
               test(array[i][j] == val_2);
            }
         }
      }

      auto size_it = array.size_begin();
      for(std::size_t i=0; i<array.size(); ++i, ++size_it) {
         test(*size_it == array[i].size());
      }
      test(size_it == array.size_end());
   } 
}

int main(int argc, char** argv)
{
   
   test_two_dimensional_variable_array<double>(1.0, 2.0);
   test_two_dimensional_variable_array<float>(1.0, 2.0);
   test_two_dimensional_variable_array<std::size_t>(1, 2);
   test_two_dimensional_variable_array<unsigned char>(1, 2); 
}

