#ifndef FEM_PRIMITIVES_H
#define FEM_PRIMITIVES_H
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

template<int dim, int N>
inline VectorXd get_shape_func_at(VectorXd x) {
    VectorXd m(N);
    if constexpr (dim == 1) {
        double x0 = x[0];
        if constexpr (N == 2) {
            m <<
                0.5 - 0.5 * x0,
                0.5 + 0.5 * x0;
        } else if constexpr (N == 3) {
            m <<
                0.5 * x0 * x0 - 0.5 * x0,
                0.5 * x0 * x0 + 0.5 * x0,
                1 - x0 * x0;
        }
    } else if constexpr (dim == 2) {
        const double x0 = x[0];
        const double x1 = x[1];
        if (N == 3) {
            m <<
                1 - x0 - x1,
                x0,
                x1;
        } else if constexpr (N == 4) {
            m <<
                0.25 * (1 - x0) * (1 - x1),
                0.25 * (1 + x0) * (1 - x1),
                0.25 * (1 + x0) * (1 + x1),
                0.25 * (1 - x0) * (1 + x1);
        }
    } else if constexpr (dim == 3) {
        const double x0 = x[0];
        const double x1 = x[1];
        const double x2 = x[2];
        if constexpr (N == 8) {
            m <<
                0.125 * (1 - x0) * (1 - x1) * (1 - x2),
                0.125 * (1 + x0) * (1 - x1) * (1 - x2),
                0.125 * (1 + x0) * (1 + x1) * (1 - x2),
                0.125 * (1 - x0) * (1 + x1) * (1 - x2),
                0.125 * (1 - x0) * (1 - x1) * (1 + x2),
                0.125 * (1 + x0) * (1 - x1) * (1 + x2),
                0.125 * (1 + x0) * (1 + x1) * (1 + x2),
                0.125 * (1 - x0) * (1 + x1) * (1 + x2);
        }
    }
    return m;
}

template<int dim, int N>
inline MatrixXd get_shape_deriv_at(VectorXd x) {
    MatrixXd m(N, dim);
    if constexpr (dim == 1) {
        double x0 = x[0];
        if constexpr (N == 2) {
            m << -0.5,
                 0.5;
        } else if constexpr (N == 3) {
            m << x0 - 0.5,
                 x0 + 0.5,
                 -x0 * 2.0;
        }
    } else if constexpr (dim == 2) {
        double x0 = x[0];
        double x1 = x[1];
        if constexpr (N == 3) {
            m <<
                -1, -1,
                1, 0,
                0, 1;
        } else if constexpr (N == 4) {
            m <<
                -0.25 * (1 - x1), -0.25 * (1 - x0),
                0.25 * (1 - x1), -0.25 * (1 + x0),
                0.25 * (1 + x1), 0.25 * (1 + x0),
                -0.25 * (1 + x1), 0.25 * (1 - x0);
        }
    } else if constexpr (dim == 3) {
        double x0 = x[0];
        double x1 = x[1];
        double x2 = x[2];
        if constexpr (N == 8) {
            m <<
                -0.125 * (1 - x1) * (1 - x2), -0.125 * (1 - x0) * (1 - x2), -0.125 * (1 - x0) * (1 - x1),
                0.125 * (1 - x1) * (1 - x2), -0.125 * (1 + x0) * (1 - x2), -0.125 * (1 + x0) * (1 - x1),
                0.125 * (1 + x1) * (1 - x2), 0.125 * (1 + x0) * (1 - x2), -0.125 * (1 + x0) * (1 + x1),
                -0.125 * (1 + x1) * (1 - x2), 0.125 * (1 - x0) * (1 - x2), -0.125 * (1 - x0) * (1 + x1),
                -0.125 * (1 - x1) * (1 + x2), -0.125 * (1 - x0) * (1 + x2), 0.125 * (1 - x0) * (1 - x1),
                0.125 * (1 - x1) * (1 + x2), -0.125 * (1 + x0) * (1 + x2), 0.125 * (1 + x0) * (1 - x1),
                0.125 * (1 + x1) * (1 + x2), 0.125 * (1 + x0) * (1 + x2), 0.125 * (1 + x0) * (1 + x1),
                -0.125 * (1 + x1) * (1 + x2), 0.125 * (1 - x0) * (1 + x2), 0.125 * (1 - x0) * (1 + x1);
        }
    }
    return m;
}

#endif //FEM_PRIMITIVES_H
