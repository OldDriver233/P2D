#ifndef FEM_PRIMITIVES_H
#define FEM_PRIMITIVES_H
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

VectorXd get_shape_func_at(size_t dim, size_t N, VectorXd x) {
    VectorXd m(N);
    if (dim == 1) {
        if(N == 2) {
            m << 0.5 - 0.5 * x(0),
                 0.5 + 0.5 * x(0);
        } else if(N == 3) {
            m << 0.5 * x(0) * x(0) - 0.5 * x(0),
                 0.5 * x(0) * x(0) + 0.5 * x(0),
                 1 - x(0) * x(0);
        }
    } else if(dim == 2) {
        if(N == 3) {
            m << x(0),
                 x(1),
                 1 - x(0) - x(1);
        } else if(N == 4) {
            m << 0.25 * (1 - x(0)) * (1 - x(1)),
                 0.25 * (1 + x(0)) * (1 - x(1)),
                 0.25 * (1 + x(0)) * (1 + x(1)),
                 0.25 * (1 - x(0)) * (1 + x(1));
        }
    }
    return m;
}

VectorXd get_shape_func_at(size_t dim, size_t N, double x) {
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

MatrixXd get_shape_deriv_at(size_t dim, size_t N, VectorXd x) {
    MatrixXd m(N, dim);
    if (dim == 1) {
        if(N == 2) {
            m << -0.5,
                  0.5;
        } else if(N == 3) {
            m <<  x(0) - 0.5,
                  x(0) + 0.5,
                 -x(0) * 2.0;
        }
    } else if(dim == 2) {
        if(N == 3) {
            m <<  1,  0,
                  0,  1,
                 -1, -1;
        } else if(N == 4) {
            m << -0.25 * (1 - x(1)), -0.25 * (1 - x(0)),
                  0.25 * (1 - x(1)), -0.25 * (1 + x(0)),
                  0.25 * (1 + x(1)),  0.25 * (1 + x(0)),
                 -0.25 * (1 + x(1)),  0.25 * (1 - x(0));
        }
    }
    return m;
}

VectorXd get_shape_deriv_at(size_t dim, size_t N, double x) {
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
