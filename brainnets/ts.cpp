#include "ts.h"

// STD
#include <iostream>
#include <map>
#include <utility>


namespace brainnets {

Eigen::MatrixXf mds(const Eigen::MatrixXf &dmtx, unsigned int dimensions)
{
    unsigned int rows = dmtx.rows();
    unsigned int cols = dmtx.cols();

    Eigen::MatrixXf identity = Eigen::MatrixXf::Identity(rows, rows);
    Eigen::MatrixXf ones = Eigen::MatrixXf::Ones(rows, rows);
    Eigen::MatrixXf centering_matrix = identity - (1.0f/rows) * ones;

    Eigen::MatrixXf input_matrix = const_cast<Eigen::MatrixXf &>(dmtx);
    Eigen::MatrixXf mtx = -0.5f * ((centering_matrix * input_matrix).cwiseProduct(input_matrix) * centering_matrix);

    Eigen::EigenSolver<Eigen::MatrixXf> eig(mtx);
    Eigen::MatrixXcf eigen_values = eig.eigenvalues();
    Eigen::MatrixXcf eigen_vectors = eig.eigenvectors();

    std::cout << std::endl;
    std::cout << eigen_values << std::endl;
    std::cout << std::endl;

    //
    // Sort the eigen vectors and values.
    //
    struct magnitude_compare {
        bool operator()(const std::complex<float>& a, const std::complex<float>& b) const {
            return std::abs(a) > std::abs(b);
        }
    };
    std::multimap<std::complex<float>, unsigned int, magnitude_compare> kvp;
    for (unsigned int i=0; i<rows; ++i) {
        std::complex<float> eigenvalue = eigen_values(i);
        kvp.insert({eigenvalue, i});
    }

    std::multimap<std::complex<float>, unsigned int>::iterator iter = kvp.begin();
    for (unsigned int i=0; i<dimensions; ++i) {
        std::complex<float> eigenvalue = iter->first;
        std::advance(iter, 1);
    }

    //
    // Extract the eigen vectors and values.
    //
    Eigen::MatrixXf extracted_eigen_vectors(rows, dimensions);
    Eigen::MatrixXf extracted_eigen_values(dimensions, 1);

    iter = kvp.begin();
    for (unsigned int i=0; i<dimensions; ++i) {
        unsigned int index = iter->second;

        extracted_eigen_vectors.col(i) = eigen_vectors.col(index).real();
        extracted_eigen_values(i) = iter->first.real();


        std::advance(iter, 1);
    }

    //
    //
    //
    Eigen::MatrixXf A(dimensions, dimensions);
    A.setZero();
    A.diagonal() = extracted_eigen_values;

    // Project the data.
    Eigen::MatrixXf X = extracted_eigen_vectors * A;

    return X;

}


Eigen::MatrixXf distance_matrix(const Eigen::MatrixXf &mtx)
{
    Eigen::MatrixXf result;

    unsigned int rows = mtx.rows();

    Eigen::MatrixXf a = mtx * mtx.transpose();
    Eigen::MatrixXf d1 = a.diagonal().replicate(1, rows);
    Eigen::MatrixXf d2 = d1.transpose() - (2.0f * a);

    result = d1 + d2;
    result = result.array().sqrt();

    return result;
}


}
