#ifndef FULL_CELL_SOLVER_H
#define FULL_CELL_SOLVER_H
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include "../shaping/primitive_type.h"
#include "stiffness/stiffness_separator.h"
#include "stiffness/stiffness_anode.h"
#include "stiffness/stiffness_cathode.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class full_cell_solver{
public:
    Primitive primitive = Primitive::Line2;
    VectorXd point_coord;
    VectorXd result;
    const int iter = 10;
    const double tolerance = 1e-12;
    int an = 6, ca = 10;
    stiffness_separator sep;
    stiffness_anode anode;
    stiffness_cathode cathode;
    int step = 0;

    explicit full_cell_solver(const VectorXd& coord): point_coord(coord) {
        sep = stiffness_separator(coord, an, ca);
        anode = stiffness_anode(coord, an, ca);
        cathode = stiffness_cathode(coord, an, ca);
    }

    void calc(Eigen::Ref<MatrixXd>);
    void apply_boundary(Eigen::Ref<MatrixXd>, Eigen::Ref<MatrixXd>, Eigen::Ref<VectorXd>);

};

#endif //FULL_CELL_SOLVER_H