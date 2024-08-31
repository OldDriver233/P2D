#ifndef FULL_CELL_SOLVER_H
#define FULL_CELL_SOLVER_H
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include "../shaping/primitive_type.h"
#include "stiffness/stiffness_separator.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class full_cell_solver{
public:
    Primitive primitive = Primitive::Line2;
    VectorXd point_coord;
    VectorXd result;
    const int iter = 10;
    const double tolerance = 1e-12;
    stiffness_separator sep;

    explicit full_cell_solver(const VectorXd& coord): point_coord(coord) {
        sep = stiffness_separator(coord);
    }

    void calc(Eigen::Ref<MatrixXd>);

};

#endif //FULL_CELL_SOLVER_H