#ifndef LPMP_LP_REPARAMETRIZATION_HXX
#define LPMP_LP_REPARAMETRIZATION_HXX

#include "hash_helper.hxx"

namespace LPMP {

   enum class lp_reparametrization_mode {Anisotropic, Anisotropic2, Uniform, Undefined};

   struct lp_reparametrization {
      lp_reparametrization() : mode(lp_reparametrization_mode::Undefined), leave_percentage(0.0) {}
      lp_reparametrization(const std::string& s)
         :
         mode(get_lp_reparametrization_mode(s)),
         leave_percentage(get_lp_reparametrization_leave_percentage(s))
      {}
      lp_reparametrization(lp_reparametrization_mode _mode, double _leave_percentage)
      : mode(_mode), leave_percentage(_leave_percentage) {}

      static double get_lp_reparametrization_leave_percentage(const std::string& s)
      {
         const auto pos = s.find(":");
         if(pos != std::string::npos) {
            return std::stod(s.substr(pos+1));
         } else {
            return 0.0;
         }
      }
      static lp_reparametrization_mode get_lp_reparametrization_mode(const std::string& s)
      {
         //feenableexcept(FE_INVALID | FE_OVERFLOW);
         const std::string uniform = "uniform";
         if(s.find("anisotropic2") == 0) {
            return lp_reparametrization_mode::Anisotropic2;
         } else if(s.find("anisotropic") == 0) {
            return lp_reparametrization_mode::Anisotropic;
         } else if(s.find("uniform") == 0) {
            return lp_reparametrization_mode::Uniform;
         } else {
            throw std::runtime_error("reparametrization mode " + s + " unknown");
         }
      } 

      bool operator==(const lp_reparametrization& o) const
      {
         return mode == o.mode && leave_percentage == o.leave_percentage; 
      }

      static std::size_t hash(const lp_reparametrization lp)
      {
         return hash::hash_combine(std::hash<double>()(lp.leave_percentage), std::hash<lp_reparametrization_mode>()(lp.mode)); 
      } 

      lp_reparametrization_mode mode;
      double leave_percentage;
   };

} // namespace LPMP

// insert hash functions for lp_reparametrization into std namespace
namespace std {

   template<>
   struct hash<typename LPMP::lp_reparametrization>
   {
      typedef typename LPMP::lp_reparametrization argument_type;
      typedef std::size_t result_type;
      result_type operator()(const argument_type& s) const
      {
         return LPMP::lp_reparametrization::hash(s);
      }
   };
} 

#endif // LPMP_LP_REPARAMETRIZATION_HXX
