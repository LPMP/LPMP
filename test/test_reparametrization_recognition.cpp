#include "test.h"
#include <string>
#include "config.hxx"

using namespace LPMP;

void test_reparametrization_recognition(const std::string s, lp_reparametrization_mode mode, const double leave_percentage)
{
   lp_reparametrization repam(s);
   test(repam.mode == mode);
   test(repam.leave_percentage == leave_percentage);
}

int main(int argc, char** argv)
{
   test_reparametrization_recognition("anisotropic", lp_reparametrization_mode::Anisotropic, 0.0);
   test_reparametrization_recognition("anisotropic:0.5", lp_reparametrization_mode::Anisotropic, 0.5);
   test_reparametrization_recognition("anisotropic:0.75", lp_reparametrization_mode::Anisotropic, 0.75);

   test_reparametrization_recognition("uniform", lp_reparametrization_mode::Uniform, 0.0);
   test_reparametrization_recognition("uniform:0.5", lp_reparametrization_mode::Uniform, 0.5);
   test_reparametrization_recognition("uniform:0.75", lp_reparametrization_mode::Uniform, 0.75);
}
