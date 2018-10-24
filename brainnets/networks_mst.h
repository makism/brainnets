/**
 *
 *
 */
#ifndef NETWORKS_MST_H
#define NETWORKS_MST_H

#include "brainnets_global.h"

// STL
#include <random>
#include <cmath>
#include <tuple>

// Eigen3
#include "eigen3/Eigen/Dense"


namespace brainnets {

namespace mst {

/**
 *
 *
 */
std::tuple<Eigen::MatrixXf, Eigen::MatrixXi> BRAINNETSSHARED_EXPORT minimum_spanning_tree(const Eigen::MatrixXf& mtx);

/**
 * A convenient function to get a 2d matrix filled with a MST's links.
 *
 */
Eigen::MatrixXi BRAINNETSSHARED_EXPORT links(const Eigen::MatrixXf& mtx);

/**
 * A convenient function to get an array filled with a MST's weights.
 *
 */
Eigen::ArrayXf BRAINNETSSHARED_EXPORT weights(const Eigen::MatrixXf& mtx);

}

}

#endif // NETWORKS_MST_H
