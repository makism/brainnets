#include "distances.h"

#include <iostream>

// Eigen
#include "eigen3/Eigen/Eigenvalues"

// Local
#include "networks.h"


namespace brainnets {

namespace distances {

float spectral_euclidean(Eigen::MatrixXf &mtx1, Eigen::MatrixXf &mtx2)
{
    if (!brainnets::are_equal_size(mtx1, mtx2)) {
        return std::nanf("");
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver1(mtx1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver2(mtx2);

    Eigen::VectorXf eigs1 = solver1.eigenvalues();
    Eigen::VectorXf eigs2 = solver2.eigenvalues();

    std::cout << eigs1 << std::endl;

    eigs1.reverseInPlace();
    eigs2.reverseInPlace();

    Eigen::VectorXf diff = (eigs1 - eigs2).array().pow(2.0);
    float sqrt_sum = std::sqrt(diff.sum());

    return sqrt_sum;
}

float spectral_k_euclidean(Eigen::MatrixXf &mtx1, Eigen::MatrixXf &mtx2, int k)
{
    if (!brainnets::are_equal_size(mtx1, mtx2)) {
        return std::nanf("");
    }

    Eigen::MatrixXf L1 = brainnets::laplacian(mtx1);
    Eigen::MatrixXf L2 = brainnets::laplacian(mtx2);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver1(L1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver2(L2);

    Eigen::VectorXf eigs1 = solver1.eigenvalues();
    Eigen::VectorXf eigs2 = solver2.eigenvalues();

    eigs1.reverseInPlace();
    eigs2.reverseInPlace();

    Eigen::ArrayXf v1 = eigs1.head(k);
    Eigen::ArrayXf v2 = eigs2.head(k);

    float num = (v1 - v2).array().pow(2.0).sum();
    float denom = 1.0;

    // Normalization factor
    float min1 = v1.array().pow(2.0).sum();
    float min2 = v2.array().pow(2.0).sum();
    if (min1 < min2) {
        denom = min1;
    } else {
        denom = min2;
    }

    float dist = std::sqrt(num / denom);

    return dist;
}

}

}
