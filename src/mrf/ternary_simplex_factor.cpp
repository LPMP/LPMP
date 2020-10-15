#include "mrf/ternary_simplex_factor.h"

namespace LPMP {

SimpleTighteningTernarySimplexFactor::SimpleTighteningTernarySimplexFactor(const std::size_t _dim1, const std::size_t _dim2, const std::size_t _dim3) 
   : msg12(_dim1, _dim2, 0.0),
   msg13(_dim1, _dim3, 0.0),
   msg23(_dim2, _dim3, 0.0)
{
//   std::fill(msg12.begin(), msg12.end(), 0.0);
//   std::fill(msg13.begin(), msg13.end(), 0.0);
//   std::fill(msg23.begin(), msg23.end(), 0.0);
}

/*
   SimpleTighteningTernarySimplexFactor(const SimpleTighteningTernarySimplexFactor& o) : 
   msg12_(o.msg12_),
   msg13_(o.msg13_),
   msg23_(o.msg23_)
   {}
   void operator=(const SimpleTighteningTernarySimplexFactor& o) {
   assert(dim1() == o.dim1() && dim2() == o.dim2() && dim3() == o.dim3());
   for(std::size_t i=0; i<dim1()*dim2(); ++i) { msg12_[i] = o.msg12_[i]; }
   for(std::size_t i=0; i<dim1()*dim3(); ++i) { msg13_[i] = o.msg13_[i]; }
   for(std::size_t i=0; i<dim2()*dim3(); ++i) { msg23_[i] = o.msg23_[i]; }
   }
 */

double SimpleTighteningTernarySimplexFactor::LowerBound() const {
   double lb = std::numeric_limits<double>::infinity();
   for(std::size_t x1=0; x1<dim1(); ++x1) {
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         for(std::size_t x3=0; x3<dim3(); ++x3) {
            lb = std::min(lb, (*this)(x1,x2,x3));
         }
      }
   }
   assert(std::isfinite(lb));
   return lb;
}

double SimpleTighteningTernarySimplexFactor::EvaluatePrimal() const
{
   if(primal_[0] >= dim1() || primal_[1] >= dim2() || primal_[2] >= dim3()) {
      return std::numeric_limits<double>::infinity();
   }
   //assert((*this)(primal_[0], primal_[1], primal_[2]) < std::numeric_limits<double>::infinity());
   return (*this)(primal_[0], primal_[1], primal_[2]);
}

void SimpleTighteningTernarySimplexFactor::MaximizePotentialAndComputePrimal()
{
   if(primal_[0] == std::numeric_limits<std::size_t>::max() && primal_[1] == std::numeric_limits<std::size_t>::max() && primal_[2] == std::numeric_limits<std::size_t>::max()) {
      double cost = std::numeric_limits<double>::infinity();
      for(std::size_t x1=0; x1<dim1(); ++x1) {
         for(std::size_t x2=0; x2<dim2(); ++x2) {
            for(std::size_t x3=0; x3<dim3(); ++x3) {
               if((*this)(x1,x2,x3) < cost) {
                  cost = (*this)(x1,x2,x3);
                  primal_ = {x1,x2,x3};
               }
            }
         }
      } 
   } else if(primal_[0] < dim1() && primal_[1] < dim2() && primal_[2] == std::numeric_limits<std::size_t>::max()) {
      double cost = std::numeric_limits<double>::infinity();
      for(std::size_t x3=0; x3<dim3(); ++x3) {
         if((*this)(primal_[0],primal_[1],x3) < cost) {
            cost = (*this)(primal_[0],primal_[1],x3);
            primal_[2] = x3;
         }
      }
   } else if(primal_[0] < dim1() && primal_[1] == std::numeric_limits<std::size_t>::max() && primal_[2] < dim3()) {
      double cost = std::numeric_limits<double>::infinity();
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         if((*this)(primal_[0],x2,primal_[2]) < cost) {
            cost = (*this)(primal_[0],x2,primal_[2]);
            primal_[1] = x2;
         }
      }
   } else if(primal_[0] == std::numeric_limits<std::size_t>::max() && primal_[1] < dim2() && primal_[2] < dim3()) {
      double cost = std::numeric_limits<double>::infinity();
      for(std::size_t x1=0; x1<dim1(); ++x1) {
         if((*this)(x1,primal_[1],primal_[2]) < cost) {
            cost = (*this)(x1,primal_[1],primal_[2]);
            primal_[0] = x1;
         }
      }
   } else if(primal_[0] < dim1() && primal_[1] == std::numeric_limits<std::size_t>::max() && primal_[2] == std::numeric_limits<std::size_t>::max()) {
      double cost = std::numeric_limits<double>::infinity();
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         for(std::size_t x3=0; x3<dim3(); ++x3) {
            if((*this)(primal_[0],x2,x3) < cost) {
               cost = (*this)(primal_[0],x2,x3);
               primal_ = {primal_[0],x2,x3};
            }
         }
      } 
   } else if(primal_[0] == std::numeric_limits<std::size_t>::max() && primal_[1] < dim2() && primal_[2] == std::numeric_limits<std::size_t>::max()) {
      double cost = std::numeric_limits<double>::infinity();
      for(std::size_t x1=0; x1<dim1(); ++x1) {
         for(std::size_t x3=0; x3<dim3(); ++x3) {
            if((*this)(x1,primal_[1],x3) < cost) {
               cost = (*this)(x1,primal_[1],x3);
               primal_ = {x1,primal_[1],x3};
            }
         }
      } 
   } else if(primal_[0] == std::numeric_limits<std::size_t>::max() && primal_[1] == std::numeric_limits<std::size_t>::max() && primal_[2] < dim3()) {
      double cost = std::numeric_limits<double>::infinity();
      for(std::size_t x1=0; x1<dim1(); ++x1) {
         for(std::size_t x2=0; x2<dim2(); ++x2) {
            if((*this)(x1,x2,primal_[2]) < cost) {
               cost = (*this)(x1,x2,primal_[2]);
               primal_ = {x1,x2,primal_[2]};
            }
         }
      } 
   } else {
      assert(false);
   }
}

double SimpleTighteningTernarySimplexFactor::operator()(const std::size_t x1, const std::size_t x2, const std::size_t x3) const 
{
   return msg12(x1,x2) + msg13(x1,x3) + msg23(x2,x3);
}

/*
   double operator[](const std::size_t x) const {
   const std::size_t x1 = x / (dim2()*dim3());
   const std::size_t x2 = ( x % (dim2()*dim3()) ) / dim3();
   const std::size_t x3 = x % dim3();
   return msg12(x1,x2) + msg13(x1,x3) + msg23(x2,x3);
   }
 */

matrix<double> SimpleTighteningTernarySimplexFactor::min_marginal12() const 
{
   matrix<double> marginal(dim1(), dim2());
   for(std::size_t x1=0; x1<dim1(); ++x1) {
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         marginal(x1,x2) = std::numeric_limits<double>::infinity();
      }
   }
   for(std::size_t x1=0; x1<dim1(); ++x1) {
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         for(std::size_t x3=0; x3<dim3(); ++x3) {
            marginal(x1,x2) = std::min( marginal(x1,x2), (*this)(x1,x2,x3));
         }
      }
   }
   return marginal;
}

matrix<double> SimpleTighteningTernarySimplexFactor::min_marginal13() const 
{
   matrix<double> marginal(dim1(), dim3()); 
   for(std::size_t x1=0; x1<dim1(); ++x1) {
      for(std::size_t x3=0; x3<dim3(); ++x3) {
         marginal(x1,x3) = std::numeric_limits<double>::infinity();
      }
   }
   for(std::size_t x1=0; x1<dim1(); ++x1) {
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         for(std::size_t x3=0; x3<dim3(); ++x3) {
            marginal(x1,x3) = std::min(marginal(x1,x3), (*this)(x1,x2,x3));
         }
      }
   }
   return marginal;
}

matrix<double> SimpleTighteningTernarySimplexFactor::min_marginal23() const 
{
   matrix<double> marginal(dim2(), dim3());
   for(std::size_t x2=0; x2<dim2(); ++x2) {
      for(std::size_t x3=0; x3<dim3(); ++x3) {
         marginal(x2,x3) = std::numeric_limits<double>::infinity();
      }
   }

   for(std::size_t x1=0; x1<dim1(); ++x1) {
      for(std::size_t x2=0; x2<dim2(); ++x2) {
         for(std::size_t x3=0; x3<dim3(); ++x3) {
            marginal(x2,x3) = std::min(marginal(x2,x3),  (*this)(x1,x2,x3));
         }
      }
   }
   return marginal;
}

std::size_t SimpleTighteningTernarySimplexFactor::dim(const std::size_t d) const 
{ 
   assert(d<3);   
   switch (d) {
      case 0 : return dim1(); break;
      case 1 : return dim2(); break;
      case 2 : return dim3(); break;
      default: return 0;
   }
}

void SimpleTighteningTernarySimplexFactor::init_primal() {
   primal_[0] = std::numeric_limits<std::size_t>::max();
   primal_[1] = std::numeric_limits<std::size_t>::max();
   primal_[2] = std::numeric_limits<std::size_t>::max();
}

} // namespace LPMP
