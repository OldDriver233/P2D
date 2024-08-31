#ifndef FEM_PRIMITIVES_H
#define FEM_PRIMITIVES_H
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

VectorXd get_shape_func_at(size_t dim, size_t N, VectorXd x);

VectorXd get_shape_func_at(size_t dim, size_t N, double x);

MatrixXd get_shape_deriv_at(size_t dim, size_t N, VectorXd x);

VectorXd get_shape_deriv_at(size_t dim, size_t N, double x);

#endif //FEM_PRIMITIVES_H
