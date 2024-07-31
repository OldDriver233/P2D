//
// Created by paradisus on 24-7-19.
//

#ifndef FEM_STIFFNESS_GENERATOR_H
#define FEM_STIFFNESS_GENERATOR_H
#include <eigen3/Eigen/Dense>
#include "../shaping/primitive_type.h"
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <vector>

using Eigen::VectorXd;
using Eigen::MatrixXd;


class stiffness_generator {
public:
    size_t dim, n;
    VectorXd points;
    std::vector<MatrixXd> cached_matrix_N, cached_matrix_dN;
    std::vector<double> cached_det_J;
    MatrixXd generated_K;
    VectorXd generated_res;
    double dt, d_ref, d, r;

    stiffness_generator() = default;
    stiffness_generator(VectorXd, Primitive, double, double, double, double);

    void generate(Eigen::Ref<MatrixXd>, Eigen::Ref<MatrixXd>);
};


#endif //FEM_STIFFNESS_GENERATOR_H
