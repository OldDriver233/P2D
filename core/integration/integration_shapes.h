#ifndef FEM_INTEGRATION_SHAPES_H
#define FEM_INTEGRATION_SHAPES_H
#include <eigen3/Eigen/Dense>
#include <cmath>

using Eigen::VectorXd;

inline VectorXd get_integration_point(size_t dim, size_t N) {
    VectorXd m(N);
    // for 1d situations
    if(dim == 1) {
        if (N == 2) {
            // line, 2 pts
            m << -sqrt(1.0 / 3.0), sqrt(1.0 / 3.0);
        } else if(N == 3) {
            // line, 3 pts
            m << -sqrt(0.6), 0, sqrt(0.6);
        }
    }
    return m;
}

inline VectorXd get_integration_weight(size_t dim, size_t N) {
    VectorXd m(N);
    if(dim == 1) {
        if (N == 2) {
            m << 1, 1;
        } else if(N == 3) {
            m << 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0;
        }
    }
    return m;
}

#endif //FEM_INTEGRATION_SHAPES_H
