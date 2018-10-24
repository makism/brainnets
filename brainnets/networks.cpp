#include "networks.h"

// STD
#include <fstream>
#include <vector>
#include <utility>
#include <map>
#include <algorithm>
#include <iterator>
#include <iostream>

// Eigen3
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include "eigen3/unsupported/Eigen/MatrixFunctions"

// Boost
#include "boost/graph/undirected_graph.hpp"
#include "boost/graph/exterior_property.hpp"
#include "boost/graph/floyd_warshall_shortest.hpp"


#include <thread>
#if !defined(WIN32)
#include <omp.h>
#include <sched.h>
#endif


typedef float t_weight;
typedef boost::property<boost::edge_weight_t, t_weight> EdgeWeightProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeWeightProperty> Graph;
typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;
typedef boost::exterior_vertex_property<Graph, t_weight> DistanceProperty;
typedef DistanceProperty::matrix_type DistanceMatrix;
typedef DistanceProperty::matrix_map_type DistanceMatrixMap;


namespace brainnets {

Eigen::MatrixXf shortest_paths(const Eigen::MatrixXf &adjacency_matrix, bool autoinverse)
{
    Eigen::MatrixXf result;

    unsigned int m = adjacency_matrix.rows();
    unsigned int n = adjacency_matrix.cols();

    //
    // Check if matrixi is square.
    //
    if (m != n) {
        return result;
    }

    //
    // Construct a graph from the adjacency matrix.
    //
    Graph graph;

    for (unsigned int x=0; x<m; ++x) {
        for (unsigned int y=x; y<n; ++y) {
            boost::add_edge(x, y, 1.0f / adjacency_matrix(x, y), graph);
        }
    }

    WeightMap weight_pmap = boost::get(boost::edge_weight, graph);
    DistanceMatrix distances(num_vertices(graph));
    DistanceMatrixMap dm(distances, graph);

    //
    // Compute all distances.
    //
    bool valid = floyd_warshall_all_pairs_shortest_paths(graph, dm, boost::weight_map(weight_pmap));

    if (!valid) {
        return result;
    }

    //
    // Prepare result matrix.
    //
    result = Eigen::MatrixXf(m, n);
    result.setZero();

    for (unsigned int x=0; x<m; ++x) {
        for (unsigned int y=x; y<n; ++y) {
            if (autoinverse) {
                result(x, y) = 1.0f / distances[x][y];
            } else {
                result(x, y) = distances[x][y];
            }
        }
    }

    result.diagonal() = Eigen::MatrixXf::Ones(m, 1);

    std::cout << result << std::endl;

    return result;
}

float global_efficiency(const Eigen::MatrixXf &mtx)
{
    unsigned int m = mtx.rows();
    unsigned int connections = (m * m - m);

    // We have to subtract m because we filled the diagonal with ones.
    float geff = (mtx.sum() - m) / static_cast<float>(connections);

    return geff;
}

Eigen::ArrayXf nodal_global_efficiency(const Eigen::MatrixXf &mtx)
{
    unsigned int m = mtx.rows();

    Eigen::MatrixXf distances = shortest_paths(mtx, true);
    Eigen::MatrixXf symm_distances = make_symmetric(distances, false);
    Eigen::ArrayXf nge = symm_distances.rowwise().sum() / (float(m - 1.0));

    return nge;
}

float cuberoot(float x)
{
    return std::pow(x, 1.0f/3);
}

Eigen::ArrayXf local_efficiency(const Eigen::MatrixXf &mtx)
{
    Eigen::ArrayXf result;

    unsigned int m = mtx.rows();
    unsigned int n = mtx.cols();

    //
    //
    //
    if (m != n) {
        return result;
    }

//    mtx.diagonal() = Eigen::MatrixXf::Ones(m, 1);

    //
    //
    //
    result = Eigen::ArrayXf(m, 1);
    result.setZero();

    Eigen::MatrixXi adjacency_matrix = adjacency_binary_matrix(mtx);

//    unsigned int nthreads = std::thread::hardware_concurrency();
//    omp_set_num_threads(nthreads - 1);

//    #pragma omp for
    for (unsigned int i=0; i<m; ++i) {
        // Find adjacent neighbors
        Eigen::MatrixXf rows = mtx.row(i);
        Eigen::MatrixXf cols = mtx.col(i);

        Eigen::MatrixXf pairs(2, m);
        pairs.row(0) = rows;
        pairs.row(1) = cols.transpose();

//        std::cout << "Pairs"
//                  << std::endl
//                  << pairs
//                  << std::endl;

        // Seperate only the pairs that at least one of two if greater than 0.0
        std::vector<int> neighbors_pairs;
        for (unsigned int ij=0; ij<m; ++ij) {
            if (pairs(0, ij) > 0.0 || pairs(1, ij) > 0.0) {
                neighbors_pairs.push_back(ij);
            }
        }

        unsigned int num_neighbors = neighbors_pairs.size();

        // Local adjacency matrix (k)
        Eigen::MatrixXi local_adjacency_matrix(num_neighbors, 1);

        // Construct a matrix containing the values of the given input
        // matrix according the found pairs
        Eigen::MatrixXf neighbors_mtx(num_neighbors, num_neighbors);

        unsigned int v1 = 0;
        unsigned int v2 = 0;
        for (unsigned int x=0; x<num_neighbors; ++x) {
            for (unsigned int y=0; y<num_neighbors; ++y) {
                v1 = neighbors_pairs[x];
                v2 = neighbors_pairs[y];
                neighbors_mtx(x, y) = mtx(v1, v2);
            }
        }
        local_adjacency_matrix = adjacency_matrix.row(i) + adjacency_matrix.col(i).transpose();

        Eigen::MatrixXf neighbors_mtx_sum(num_neighbors, 1);
        neighbors_mtx_sum.setZero();
        unsigned int offset = 0;
        for (std::vector<int>::const_iterator it = neighbors_pairs.begin();
             it != neighbors_pairs.end();
             ++it, ++offset)
        {
            neighbors_mtx_sum(offset, 0) = std::pow(pairs(0, *it), 1.0f/3) + std::pow(pairs(1, *it), 1.0f/3);
        }

        // Shortest paths among the neighbors
        Eigen::MatrixXf neighbors_di = shortest_paths(neighbors_mtx).unaryExpr(&cuberoot);
        neighbors_di.triangularView<Eigen::Lower>() = neighbors_di.transpose();
        neighbors_di.diagonal() = Eigen::MatrixXf::Zero(num_neighbors, 1);

        // Note:
        // Convert to .array() to do element-wise multiplication
        Eigen::MatrixXf weights_paths = ((neighbors_mtx_sum * neighbors_mtx_sum.transpose())).array() * neighbors_di.array();
        float coeff = weights_paths.sum();

        // We have to subtract 2, because of the adjacent matrix, the self loops indicate 2 because of the symmetrization
        float degrees_coeff = std::pow(local_adjacency_matrix.sum(), 2.0f) - local_adjacency_matrix.sum() - 2;

        std::cout << "Numerator \t" << coeff << std::endl;
        std::cout << "Denumerator \t" << degrees_coeff << std::endl;

        result[i] = 0.5f * coeff / degrees_coeff;
    }


    return result;
}

Eigen::ArrayXf degree(const Eigen::MatrixXf &mtx)
{
    Eigen::MatrixXi binary = (mtx.array() != 0.0f).cast<int>();
    Eigen::ArrayXf result = binary.rowwise().sum().cast<float>();

    return result;
}

Eigen::MatrixXi adjacency_binary_matrix(const Eigen::MatrixXf &mtx, bool symmetrize)
{
    unsigned int m = mtx.rows();

    Eigen::MatrixXi adjacency_matrix = (mtx.array() > 0.0f).cast<int>();
    if (symmetrize) {
        adjacency_matrix.triangularView<Eigen::Lower>() = adjacency_matrix.transpose();
    }
    adjacency_matrix.diagonal() = Eigen::MatrixXi::Ones(m, 1);

    return adjacency_matrix;
}

Eigen::MatrixXf watts_strogatz(unsigned int n, unsigned int k, float p)
{
    // Prepare random generator
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<unsigned int> n_dist{0, n - 1};
    std::uniform_real_distribution<float> p_dist{0.0, 1.0};

    // Prepare matrices
    Eigen::MatrixXi adj = Eigen::MatrixXi(n, n);
    adj.setZero();

    Eigen::MatrixXf mtx = Eigen::MatrixXf(n, n);
    mtx.setZero();

    Eigen::ArrayXi sum_per_row = Eigen::ArrayXi(n);
    sum_per_row.setConstant(k + 1);

    // Fill in the initial neighbors
    for (unsigned int x=0; x<n; ++x) {
        adj(x, x) = 0;

        std::vector<int> neighbors;

        int x1 = -(k / 2.0);
        int x2 = k / 2.0;
        for (int neighbor = x1; neighbor <= x2; ++neighbor) {
            if (neighbor == 0)
                continue;

            int offset = (int)x + neighbor;


            if (offset < 0)
                offset += n;
            if (offset >= n)
                offset -= n;

            neighbors.push_back(offset);
        }

        for (unsigned int neighbor=0; neighbor<neighbors.size(); ++neighbor) {
            unsigned int nn = neighbors[neighbor];
            adj(x, nn) = 1;
        }
    }

    // Apply model
    for (unsigned int i=1; i<k-1; ++i) {
        for (unsigned int node=0; node<n; ++node) {
            int n1 = node + i;
            int n2 = node - i;

            if (n1 < 0) n1 += n;
            if (n1 >= n) n1 -= n;

            if (n2 < 0) n2 += n;
            if (n2 >= n) n2 -= n;


            // Neighbor 11
            float r = p_dist(generator);
            if (r < p) {
                unsigned int w = n_dist(generator);

                if (adj(node, w) == 0) {
                    adj(node, w) = 1;
                    adj(node, n1) = 0;
                }
            }

            // Neighbor #2
            r = p_dist(generator);
            if (r < p) {
                unsigned int w = n_dist(generator);

                if (adj(node, w) == 0) {
                    adj(node, w) = 1;
                    adj(node, n2) = 0;
                }
            }


        }
    }

    Eigen::ArrayXi curr_row_per_sum = adj.rowwise().sum();
    if ((sum_per_row == curr_row_per_sum).all()) {

    }


    mtx = adj.cast<float>();

    return mtx;
}

Eigen::ArrayXf clustering_coeff(const Eigen::MatrixXf &mtx)
{
    unsigned int m = mtx.rows();
    unsigned int n = mtx.cols();

    Eigen::ArrayXf k = degree(mtx);

    Eigen::MatrixXf result = mtx;
    for (unsigned int x=0; x<m; ++x) {
        for (unsigned int y=0; y<n; ++y) {
            result(x, y) = std::pow(result(x, y), 1.0/3.0);
        }
    }

    result = result * result * result;

    Eigen::ArrayXf diag = result.diagonal();
    Eigen::ArrayXf result2 = diag.cwiseQuotient(k.cwiseProduct(k - 1));

    return result2;
}

bool is_symmetric(const Eigen::MatrixXf &matrix)
{
    Eigen::MatrixXf upper = matrix.triangularView<Eigen::StrictlyUpper>();
    Eigen::MatrixXf lower = matrix.triangularView<Eigen::StrictlyLower>();
    float sum = (upper - lower).sum();

    if (sum != 0.0f) {
        return false;

    } else {
        return true;

    }
}

Eigen::MatrixXf make_symmetric(const Eigen::MatrixXf &matrix, bool fillDiagonalWithOnes)
{
    Eigen::MatrixXf result(matrix);
    unsigned int m = matrix.rows();

    result.triangularView<Eigen::Lower>() = matrix.transpose();
    if (fillDiagonalWithOnes) {
        result.diagonal() = Eigen::MatrixXf::Ones(m, 1);
    } else {
        result.diagonal() = Eigen::MatrixXf::Zero(m, 1);
    }

    return result;

}

Eigen::MatrixXf triangles(const Eigen::MatrixXf &mtx)
{
    Eigen::MatrixXf result;

    return result;

}

Eigen::ArrayXf strength(const Eigen::MatrixXf &mtx)
{
    return mtx.colwise().sum();
}

Eigen::ArrayXf closeness_centrality(const Eigen::MatrixXf &mtx)
{
    Eigen::ArrayXf result;


    return result;
}

Eigen::ArrayXf within_module_degree_z_score(const Eigen::MatrixXf &mtx)
{
    Eigen::ArrayXf result;

    return result;
}

Eigen::MatrixXf filter_matrix(const Eigen::MatrixXf &mtx, float min, float max)
{
    Eigen::MatrixXf pass1 = (mtx.array() <= max).select(mtx, 0.0f);
    Eigen::MatrixXf pass2 = (pass1.array() >= min).select(pass1, 0.0f);

    return pass2;
}

Eigen::MatrixXf filter_percentage_matrix(const Eigen::MatrixXf &mtx, float percetange, bool skip_diagonal)
{
    int rows = mtx.rows();
    int cols = mtx.cols();

    if (rows != cols) {
        return mtx;
    }

    Eigen::MatrixXf result(rows, cols);
    result.setZero();

    std::multimap<float, int, std::greater<float>> kvp;
    std::map<int, std::pair<int, int>> indices_to_coords;

    int indice = 0;
    for (int x=0; x<rows; ++x) {
        for (int y=0; y<cols; ++y) {
            if (x == y && skip_diagonal)
                continue;

            kvp.insert({mtx(x, y),indice});
            indices_to_coords[indice] = std::make_pair(x, y);

            indice++;
        }
    }

    int total_elements = rows * cols;
    int elements = std::floor(total_elements * percetange);

    std::multimap<float, int>::iterator it = kvp.begin();
    for (int i=0; i<elements; ++i) {
        std::pair<int, int> xy = indices_to_coords[it->second];

        result(xy.first, xy.second) = it->first;

        std::advance(it, 1);
    }

    return result;
}

bool to_csv(const Eigen::MatrixXf &matrix, const std::string &name)
{
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream file(name.c_str());
    file << matrix.format(CSVFormat);

    file.close();

    return true;
}

void treatment(Eigen::MatrixXf &mtx)
{
    int rows = mtx.rows();
    int cols = mtx.cols();

    for (int x=0; x<rows; ++x) {
        for (int y=0; y<cols; ++y) {
            float val = mtx(x, y);

            if (std::isnan(val) || std::isinf(val)) {
                mtx(x, y) = 0.0f;
            }
        }
    }

}

Eigen::ArrayXf eigenvector_centrality(const Eigen::MatrixXf &mtx)
{
    Eigen::EigenSolver<Eigen::MatrixXf> es(mtx, true);
    Eigen::MatrixXcf eigv = es.eigenvectors();

    std::cout << std::endl;
    std::cout << es.eigenvalues() << std::endl;
    std::cout << std::endl;

    return eigv.col(0).real().array().abs();
}

Eigen::MatrixXf laplacian(const Eigen::MatrixXf &mtx, bool normed)
{
    Eigen::MatrixXf L = Eigen::MatrixXf::Zero(mtx.rows(), mtx.cols());
    Eigen::MatrixXf D = brainnets::degree_weighted_undirected(mtx);
    L.diagonal() <<  D;
    L = L - mtx;

    return L;
}

bool are_equal_size(const Eigen::MatrixXf &mtx1, const Eigen::MatrixXf &mtx2)
{
    int rows1 = mtx1.rows();
    int cols1 = mtx1.cols();

    int rows2 = mtx2.rows();
    int cols2 = mtx2.cols();

    return (rows1 == rows2 && cols1 == cols2);
}

Eigen::ArrayXf degree_weighted_undirected(const Eigen::MatrixXf &mtx)
{
    return mtx.rowwise().sum();
}

}
