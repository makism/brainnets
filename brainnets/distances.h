#ifndef DISTANCES_H
#define DISTANCES_H

#include "brainnets_global.h"

// STL

// Eigen3
#include "eigen3/Eigen/Dense"

namespace brainnets {

namespace distances {

/**
 * @brief spectral_euclidean
 * @param mtx1
 * @param mtx2
 * @return
 */
float BRAINNETSSHARED_EXPORT spectral_euclidean(Eigen::MatrixXf &mtx1, Eigen::MatrixXf &mtx2);

/**
 * @brief spectral_k_euclidean
 * @param mtx1
 * @param mtx2
 * @param k
 * @return
 */
float BRAINNETSSHARED_EXPORT spectral_k_euclidean(Eigen::MatrixXf &mtx1, Eigen::MatrixXf &mtx2, int k=3);

}

}

#endif // DISTANCES_H
