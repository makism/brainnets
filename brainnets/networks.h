/**
 * Networks
 *
 *
 *
 */
#ifndef NETWORKS_H
#define NETWORKS_H

#include "brainnets_global.h"

// STL
#include <random>
#include <cmath>

// Eigen3
#include "eigen3/Eigen/Dense"


namespace brainnets {

/**
 * Equal size.
 *
 * Check if two matrices have the same dimensions.
 *
 *
 * @param   mtx1  Matrix 1
 * @param   mtx2  Matrix 2
 * @return  boolean
 */
bool BRAINNETSSHARED_EXPORT are_equal_size(const Eigen::MatrixXf &mtx1, const Eigen::MatrixXf &mtx2);


/**
 * Treament matrix.
 *
 * Replace NaN, -Inf, +Inf with 0.0f.
 *
 * The operation takes inplace.
 *
 */
void BRAINNETSSHARED_EXPORT treatment(Eigen::MatrixXf& mtx);


/**
 *
 *
 *
 */
bool BRAINNETSSHARED_EXPORT to_csv(const Eigen::MatrixXf& matrix, const std::string& name);


/**
 *
 *
 */
bool BRAINNETSSHARED_EXPORT is_symmetric(const Eigen::MatrixXf& matrix);


/**
 *
 *
 */
Eigen::MatrixXf BRAINNETSSHARED_EXPORT make_symmetric(const Eigen::MatrixXf& matrix, bool fillDiagonalWithOnes=true);


/**
 * Generate a random 'small' network based on the Watts and Strogatz method.
 *
 * @param n Number of nodes
 * @param k Number of neighbors
 * @param p Probability of rewiring a edge
 *
 */
Eigen::MatrixXf BRAINNETSSHARED_EXPORT watts_strogatz(unsigned int n, unsigned int k, float p = 0.10f);


/**
 * Clustering Coefficient
 *
 * @see https://arxiv.org/pdf/physics/0612169.pdf
 * @see https://pdfs.semanticscholar.org/da00/6d704311c4183d4b528331bac84b9bd312f5.pdf
 * @see http://www.bioinformatics.org/aminonet/doc/formu.pdf
 * @see https://toreopsahl.com/tnet/weighted-networks/clustering/
 * @see http://www.pnas.org/content/101/11/3747.full
 * @see http://jponnela.com/web_documents/a9.pdf
 * @see https://arxiv.org/pdf/cond-mat/0311416.pdf
 * @see http://journals.aps.org.secure.sci-hub.cc/pre/abstract/10.1103/PhysRevE.71.065103
 */
Eigen::ArrayXf BRAINNETSSHARED_EXPORT clustering_coeff(const Eigen::MatrixXf& mtx);


/**
 *
 *
 *
 */
Eigen::MatrixXf BRAINNETSSHARED_EXPORT laplacian(const Eigen::MatrixXf& mtx, bool normed=true);

/**
 * @brief shortest_paths
 *
 * @param       adjacency_matrix
 * @param       autoinverse         Compute the inverse distances.
 * @return
 */
Eigen::MatrixXf BRAINNETSSHARED_EXPORT shortest_paths(const Eigen::MatrixXf& adjacency_matrix, bool autoinverse=true);


/**
 * Global efficiency
 *
 * @param       mtx     Connectivity matrix
 * @return              Global efficiency (scalar)
 */
float BRAINNETSSHARED_EXPORT global_efficiency(const Eigen::MatrixXf& mtx);


/**
 *
 *
 *
 */
Eigen::ArrayXf BRAINNETSSHARED_EXPORT nodal_global_efficiency(const Eigen::MatrixXf &mtx);

/**
 * Local efficiency
 *
 * @param       mtx     Connectivity matrix
 * @return              Local efficiency (matrix)
 */
Eigen::ArrayXf BRAINNETSSHARED_EXPORT local_efficiency(const Eigen::MatrixXf& mtx);


/**
 * Degree
 *
 * @param mtx
 * @return
 */
Eigen::ArrayXf BRAINNETSSHARED_EXPORT degree(const Eigen::MatrixXf& mtx);
Eigen::ArrayXf BRAINNETSSHARED_EXPORT degree_weighted_undirected(const Eigen::MatrixXf& mtx);

/**
 * Strength
 *
 */
Eigen::ArrayXf BRAINNETSSHARED_EXPORT strength(const Eigen::MatrixXf& mtx);


/**
 * Triangles
 *
 *
 */
Eigen::MatrixXf BRAINNETSSHARED_EXPORT triangles(const Eigen::MatrixXf& mtx);


/**
 * Closeness Centrality
 *
 */
Eigen::ArrayXf BRAINNETSSHARED_EXPORT closeness_centrality(const Eigen::MatrixXf& mtx);


/**
 *
 *
 */
Eigen::ArrayXf BRAINNETSSHARED_EXPORT eigenvector_centrality(const Eigen::MatrixXf& mtx);

/**
 * Characteristic Path Length
 *
 *
 */
Eigen::ArrayXf BRAINNETSSHARED_EXPORT characteristic_path_length(const Eigen::MatrixXf& mtx);


/**
 * Within Module Degree Z-Score
 *
 *
 */
Eigen::ArrayXf BRAINNETSSHARED_EXPORT within_module_degree_z_score(const Eigen::MatrixXf& mtx);


/**
 * Adjacency binary matrix
 *
 * Given a connectivity matrix, the binary adjacency matrix is computed.
 * This is easily done:
 *  1) First, we filter out all negative values.
 *  2) Then, we just put a "one" where a value exists, and "zero" where zero.
 *  3) Finally, make the matrix symmetric.
 *
 * @aparam      mtx  Connectivity matrix
 * @return           Adjacency binary matrix
 */
Eigen::MatrixXi adjacency_binary_matrix(const Eigen::MatrixXf& mtx, bool symmetrize=false);

/**
 *
 *
 *
 */
Eigen::MatrixXf filter_matrix(const Eigen::MatrixXf& mtx, float min, float max = 1.0f);

/**
 *
 *
 *
 */
Eigen::MatrixXf filter_percentage_matrix(const Eigen::MatrixXf& mtx, float percetange = 0.30f, bool skip_diagonal = true);

}

#endif // NETWORKS_H
