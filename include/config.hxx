#pragma once

#include <stdexcept>
#include <string>
#include <array>
#include <cmath>
#include <cassert>
#include <limits>
#include "tclap/CmdLine.h"
#include "lp_reparametrization.hxx"
#include "hash_helper.hxx"

#define SIMDPP_ARCH_X86_AVX2
#include "simdpp/simd.h"

// type definitions for LPMP

namespace LPMP {

   // data types for all floating point/integer operations 
   // float is inaccurate for large problems and I observed oscillation. Possibly, some sort of numerical stabilization needs to be employed
   //using REAL = float;
   //constexpr std::size_t REAL_ALIGNMENT = 8;
   //using REAL_VECTOR = simdpp::float32<REAL_ALIGNMENT>;

   using REAL = double;
   constexpr std::size_t REAL_ALIGNMENT = 4;
   using REAL_VECTOR = simdpp::float64<REAL_ALIGNMENT>;

   using INDEX = std::size_t;
   using UNSIGNED_INDEX = INDEX;
   using SIGNED_INDEX = long int; // note: must be the same as flow type in MinCost
   using SHORT_INDEX = unsigned char;
   using LONG_SIGNED_INDEX = long int;
   using LONG_INDEX = long unsigned int;

   enum class Chirality {left,right};
   enum class inequality_type { smaller_equal, greater_equal, equal };

   enum class MessageSendingType {SRMP,MPLP}; // TODO: remove
   enum class Direction {forward, backward};

   constexpr REAL eps = std::is_same<REAL,float>::value ? 1e-6 : 1e-8;
   // verbosity levels: 0: silent
   //                   1: important diagnostics, e.g. lower bound, upper bound, runtimes
   //                   2: debug informations
   static std::size_t verbosity = 0; 
   static bool diagnostics() { return verbosity >= 1; }
   static bool debug() { return verbosity >= 2; }
   
   // steers optimization of LP solver. Is returned by visitor and processed by solver.
   // also put this into solver.hxx
   class LpControl {
   public:
      lp_reparametrization lp_repam;
      bool computePrimal = false;
      bool computeLowerBound = false;
      bool tighten = false;
      bool end = false; // terminate optimization
      bool error = false;
      std::size_t tightenConstraints = 0; // when given as return type, indicates how many constraints are to be added. When given as parameter to visitor, indicates how many were added.
      double tightenMinDualIncrease = 0.0; // do zrobienia: obsolete
   }; 

   template<typename T>
   T normalize(const T x) {
      static_assert(std::is_same<T,double>::value || std::is_same<T,float>::value,"");
      assert(!std::isnan(x));
      if(std::isfinite(x)) {
         return x;
      } else {
         return std::numeric_limits<REAL>::infinity();
      }
   }

   // TCLAP constraints
   class PositiveRealConstraint : public TCLAP::Constraint<REAL>
   {
      public:
         std::string description() const { return "positive real constraint"; };
         std::string shortID() const { return "positive real number"; };
         bool check(const REAL& value) const { return value >= 0.0; };
   };
   class OpenUnitIntervalConstraint: public TCLAP::Constraint<REAL>
   {
      public:
         std::string description() const { return "0<x<1 real constraint"; };
         std::string shortID() const { return "positive real number smaller 1"; };
         bool check(const REAL& value) const { return value > 0.0 && value < 1.0; };
   };
   class ClosedUnitIntervalConstraint: public TCLAP::Constraint<REAL>
   {
      public:
         std::string description() const { return "0<=x<=1 real constraint"; };
         std::string shortID() const { return "non-negative real number smaller equal 1"; };
         bool check(const REAL& value) const { return value >= 0.0 && value <= 1.0; };
   };
   class PositiveIntegerConstraint : public TCLAP::Constraint<INDEX>
   {
      public:
         std::string description() const { return "strictly positive integer constraint"; };
         std::string shortID() const { return "strictly positive integer"; };
         bool check(const INDEX& value) const { return value > 0; };
   };
   static PositiveIntegerConstraint positiveIntegerConstraint;

} // namespace LPMP
