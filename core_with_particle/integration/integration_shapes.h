#ifndef FEM_INTEGRATION_SHAPES_H
#define FEM_INTEGRATION_SHAPES_H
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <numbers>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::numbers::inv_sqrt3_v;

template<int dim, int N>
inline int get_integration_point_count() {
    int ret;
    ret = N;
    return ret;
}

template<int dim, int N>
inline MatrixXd get_integration_point() {
    int pts = get_integration_point_count<dim, N>();
    MatrixXd m(pts, dim);
    // for 1d situations
    if constexpr (dim == 1) {
        if constexpr (N == 2) {
            // line, 2 pts
            m <<
                -inv_sqrt3_v<double>,
                inv_sqrt3_v<double>;
        } else if constexpr (N == 3) {
            // line, 3 pts
            m <<
                -sqrt(0.6),
                0,
                sqrt(0.6);
        }
    } else if constexpr (dim == 2) {
        if constexpr (N == 3) {
            // 3pt Hammer integration
            m <<
                0.2, 0.2,
                0.6, 0.2,
                0.2, 0.6;
        } else if constexpr (N == 4) {
            m <<
                -inv_sqrt3_v<double>, -inv_sqrt3_v<double>,
                inv_sqrt3_v<double>, -inv_sqrt3_v<double>,
                inv_sqrt3_v<double>, inv_sqrt3_v<double>,
                -inv_sqrt3_v<double>, inv_sqrt3_v<double>;
        }
    }
    return m;
}

template<int dim, int N>
inline VectorXd get_integration_weight() {
    int pts = get_integration_point_count<dim, N>();
    VectorXd m(pts);
    if constexpr (dim == 1) {
        if constexpr (N == 2) {
            m <<
                1,
                1;
        } else if constexpr (N == 3) {
            m <<
                5. / 9.,
                8. / 9.,
                5. / 9.;
        }
    } else if constexpr (dim == 2) {
        if constexpr (N == 3) {
            m <<
                1. / 6.,
                1. / 6.,
                1. / 6.;
        } else if constexpr (N == 4) {
            m <<
                1,
                1,
                1,
                1;
        }
    }
    return m;
}

#endif //FEM_INTEGRATION_SHAPES_H
