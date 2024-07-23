//
// Created by paradisus on 24-7-19.
//

#ifndef FEM_PARTICLE_SOLVER_H
#define FEM_PARTICLE_SOLVER_H
#include <eigen3/Eigen/Dense>
#include <utility>
#include "../shaping/primitive_type.h"
#include "stiffness_generator.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;


class particle_solver {
public:
    Primitive primitive = Primitive::Line2;
    VectorXd point_coord;
    VectorXd result;
    stiffness_generator gen;
    const int iter = 10;
    const double tolerance = 1e-12;
    const double dt = 1.0;
    const double d_ref = 3.9e-14;
    const double d = 3.9e-14;
    const double r = 1e-5;

    explicit particle_solver(const VectorXd& coord): point_coord(coord) {
        gen = stiffness_generator(coord, primitive, dt, d_ref, d, r);
    }

    void calc(Eigen::Ref<MatrixXd>);

};


#endif //FEM_PARTICLE_SOLVER_H
