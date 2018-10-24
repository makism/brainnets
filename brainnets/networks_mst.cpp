#include "networks_mst.h"

// STD
#include <iostream>
#include <vector>
#include <tuple>

// Boost
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/prim_minimum_spanning_tree.hpp"
#include "boost/graph/kruskal_min_spanning_tree.hpp"

// Local
#include "ts.h"


namespace brainnets {

namespace mst {

std::tuple<Eigen::MatrixXf, Eigen::MatrixXi> minimum_spanning_tree(const Eigen::MatrixXf &mtx)
{
    unsigned int rows = mtx.rows();
    unsigned int cols = mtx.cols();

    //
    //
    //
    Eigen::MatrixXf dmtx = distance_matrix(mtx);

    //
    //
    //
    typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,
            boost::no_property, boost::property < boost::edge_weight_t, float > > Graph;
    typedef boost::graph_traits < Graph >::edge_descriptor Edge;
    typedef std::pair<int, int> E;

#if defined(WIN32)
    const int num_nodes = rows * cols;
    E *edge_array = new E[num_nodes]();
    float *weights = new float[num_nodes]();

    int counter = 0;
    for (int i=0; i<rows; ++i) {
        for (int ii=0; ii<cols; ++ii) {
            edge_array[counter] = E(i, ii);
            weights[counter++] = dmtx(i, ii);
        }
    }

    std::size_t num_edges = rows * cols;
    std::cout << "num_edges =" << num_edges << std::endl;
#else
    const int num_nodes = rows * cols;
    E edge_array[rows * cols] = {};
    float weights[rows * cols] = {};

    int counter = 0;
    for (int i=0; i<rows; ++i) {
        for (int ii=0; ii<cols; ++ii) {
            edge_array[counter] = E(i, ii);
            weights[counter++] = dmtx(i, ii);
        }
    }

    std::size_t num_edges = sizeof(edge_array) / sizeof(E);
#endif

    Graph g(edge_array, edge_array + num_edges, weights, num_nodes);


    boost::property_map < Graph, boost::edge_weight_t >::type weight = boost::get(boost::edge_weight, g);
    std::vector<Edge> spanning_tree;

    boost::kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

    //
    //
    //
    Eigen::MatrixXf mst(rows, cols);
    mst.setZero();

    Eigen::MatrixXi bin_mst(rows, cols);
    bin_mst.setZero();

    for (std::vector<Edge>::iterator iter = spanning_tree.begin();
         iter != spanning_tree.end();
         ++iter)
    {
        int s = source(*iter, g);
        int e = target(*iter, g);
        float w = weight[*iter];

        mst(s, e) = w;
        bin_mst(s, e) = 1;
    }

    //
    // Prepare result
    //
    std::tuple<Eigen::MatrixXf, Eigen::MatrixXi> mst_results =
            std::make_tuple(mst, bin_mst);

    return mst_results;
}

Eigen::MatrixXi links(const Eigen::MatrixXf &mtx)
{
    unsigned int num_links = (mtx.array() != 0.0f).count();
    unsigned int rows = mtx.rows();
    unsigned int cols = mtx.cols();

    Eigen::MatrixXi result(num_links, 2);

    unsigned int cntr = 0;
    for (int i=0; i<rows; ++i) {
        for (int ii=0; ii<cols; ++ii) {
            float val = mtx(i, ii);

            if (val != 0.0f) {
                result(cntr, 0) = i;
                result(cntr, 1) = ii;

                cntr++;
            }
        }
    }

    return result;
}

Eigen::ArrayXf weights(const Eigen::MatrixXf &mtx)
{
    unsigned int num_links = (mtx.array() != 0.0f).count();
    unsigned int rows = mtx.rows();
    unsigned int cols = mtx.cols();

    Eigen::ArrayXf result(num_links);

    unsigned int cntr = 0;
    for (int i=0; i<rows; ++i) {
        for (int ii=0; ii<cols; ++ii) {
            float val = mtx(i, ii);

            if (val != 0.0f) {
                result(cntr) = val;

                cntr++;
            }
        }
    }

    return result;
}

}

}
