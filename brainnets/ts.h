/**
 * Timeseries
 *
 * @author Avraam Marimpis <marimpis@brainvoyager.com
 *
 */

#ifndef TS_H
#define TS_H

#include "brainnets_global.h"

// Eigen3
#include "eigen3/Eigen/Dense"


namespace brainnets {

void surrogates();

void fdr();

void aaft();

void phase_rand();


/**
 * Compute a distance matrix using the Euclidean distance.
 *
 * rows = recordings
 * cols = features (time?)
 *
 */
Eigen::MatrixXf distance_matrix(const Eigen::MatrixXf& mtx);


/**
 * Classical (aka Metric) Multidimensional Scaling.
 *
 *
 */
Eigen::MatrixXf mds(const Eigen::MatrixXf& dmtx, unsigned int dimensions = 2);

}

#endif // TS_H
