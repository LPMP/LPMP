#include <Eigen/Eigen>

namespace LPMP {

    Eigen::SparseMatrix<double> commutation_matrix(const std::size_t n, const std::size_t m)
    {
        Eigen::SparseMatrix<double> K(n*m,n*m);
        K.reserve(Eigen::VectorXi::Constant(K.cols(), 1));

        for(std::size_t i=0; i<n; ++i) {
            for(std::size_t j=0; j<m; ++j) {
                K.insert(j*n+i, i*m + j) = 1.0;
            }
        }

        return K;
    }

    Eigen::SparseMatrix<double> identity_matrix(const std::size_t n)
    {
        Eigen::SparseMatrix<double> I(n,n);
        I.setIdentity();
        return I; 
    }
}
