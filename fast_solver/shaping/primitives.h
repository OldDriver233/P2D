#ifndef FEM_PRIMITIVES_H
#define FEM_PRIMITIVES_H
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

template<int dim, int N>
inline VectorXd get_shape_func_at(double x) {
    VectorXd m(N);
    if (dim == 1) {
        if(N == 2) {
            m << 0.5 - 0.5 * x,
                 0.5 + 0.5 * x;
        } else if(N == 3) {
            m << 0.5 * x * x - 0.5 * x,
                 0.5 * x * x + 0.5 * x,
                 1 - x * x;
        }
    }
    return m;
}

template<int dim, int N>
inline VectorXd get_shape_deriv_at(double x) {
    VectorXd m(N);
    if (dim == 1) {
        if(N == 2) {
            m << -0.5,
                 0.5;
        } else if(N == 3) {
            m << x - 0.5,
                 x + 0.5,
                 -x * 2.0;
        }
    }
    return m;
}

#endif //FEM_PRIMITIVES_H
