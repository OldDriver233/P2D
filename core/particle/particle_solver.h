//
// Created by paradisus on 24-7-19.
//

#ifndef FEM_PARTICLE_SOLVER_H
#define FEM_PARTICLE_SOLVER_H
#include <eigen3/Eigen/Dense>
#include <utility>
#include "../shaping/primitive_type.h"
#include "stiffness_generator.h"
#include "../constants/constant.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;


class particle_solver {
public:
    Primitive primitive = Primitive::Line2;
    VectorXd point_coord;
    VectorXd result;
    stiffness_generator gen;
    const int iter = 10;
    const double tolerance = constant::tolerance;
    const double dt = constant::dt;
    const double r = constant::r;
    int type = constant::type;

    explicit particle_solver(const VectorXd& coord): point_coord(coord) {
        gen = stiffness_generator(coord, primitive, dt, r);
    }

    void calc(Eigen::Ref<MatrixXd>);

};


#endif //FEM_PARTICLE_SOLVER_H
