#include "commutation_matrix.h"
#include <eigen3/Eigen/Eigen>
#include "test.h"

using namespace LPMP;

int main(int argc, char** argv)
{

    for(std::size_t i=1; i<20; ++i) {
        for(std::size_t j=1; j<20; ++j) {
            const auto K = commutation_matrix(i,j);
            const auto Id = Eigen::MatrixXd::Identity(i*j,i*j);
            test( (K*K.transpose() - Id).norm() <= 1e-8);
            test( (K.transpose()*K - Id).norm() <= 1e-8);

            // get random matrix and check whether K permutes
            const Eigen::MatrixXd M = Eigen::MatrixXd::Random(i,j);
            const Eigen::Map<const Eigen::VectorXd> M_v(M.data(), M.size());
            const Eigen::MatrixXd Mt = M.transpose();
            const Eigen::Map<const Eigen::VectorXd> Mt_v(Mt.data(), Mt.size());

            test( (M_v - K*Mt_v).norm() <= 1e-8);

        }
    }

}
