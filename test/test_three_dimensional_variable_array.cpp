#include "test.h"
#include "three_dimensional_variable_array.hxx"
#include <random>

using namespace LPMP;

template<typename T>
void test_three_dimensional_variable_array(T val_1, T val_2)
{
   // random vector tests
   std::random_device rd{};

   for(std::size_t n=2; n<100; ++n) {
      std::mt19937 gen{rd()};
      std::uniform_int_distribution<> dist{0,40};
      std::vector<std::array<std::size_t,2>> size(n);
      for(auto& x : size) {
         x[0] = dist(gen);
         x[1] = dist(gen);
      }

      three_dimensional_variable_array<T> array(size.begin(), size.begin() + size.size()/2, val_1);
      test(array.size() == size.size()/2);
      for(std::size_t i=0; i<array.size(); ++i) {
         test(array[i].size() == size[i][0]*size[i][1]); 
         test(array[i].dim1() == size[i][0]);
         test(array[i].dim2() == size[i][1]);
      }
      for(std::size_t i=0; i<array.size(); ++i) { 
         for(std::size_t j=0; j<array[i].dim1(); ++j) { 
            for(std::size_t k=0; k<array[i].dim2(); ++k) { 
               test(array(i,j,k) == val_1);
            }
         }
      }

      array.resize(size.begin(), size.end(), 2.0);
      test(array.size() == size.size());
      for(std::size_t i=0; i<array.size(); ++i) { test(array[i].size() == size[i][0]*size[i][1]); }
      for(std::size_t i=0; i<array.size(); ++i) { 
         for(std::size_t j=0; j<array[i].dim1(); ++j) { 
            for(std::size_t k=0; k<array[i].dim2(); ++k) { 
               if(i < size.size()/2) {
                  test(array(i,j,k) == val_1);
               } else {
                  test(array(i,j,k) == val_2);
               }
            }
         }
      }

      auto size_it = array.size_begin();
      for(std::size_t i=0; i<array.size(); ++i, ++size_it) {
         test((*size_it)[0] == array[i].dim1());
         test((*size_it)[1] == array[i].dim2());
      }
      test(size_it == array.size_end()); 
   } 
}

int main(int argc, char** argv)
{
   
   test_three_dimensional_variable_array<double>(1.0, 2.0);
   test_three_dimensional_variable_array<float>(1.0, 2.0);
   test_three_dimensional_variable_array<std::size_t>(1, 2);
   test_three_dimensional_variable_array<unsigned char>(1, 2); 
}
