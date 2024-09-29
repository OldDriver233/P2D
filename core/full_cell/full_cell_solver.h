#ifndef FULL_CELL_SOLVER_H
#define FULL_CELL_SOLVER_H
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseLU>
//#include <eigen3/Eigen/PardisoSupport>
//#include <eigen3/Eigen/UmfPackSupport>
//#include <eigen3/Eigen/src/UmfPackSupport/UmfPackSupport.h>
#include "../shaping/primitive_type.h"
#include "stiffness/stiffness_separator.h"
#include "stiffness/stiffness_anode.h"
#include "stiffness/stiffness_cathode.h"
#include "../constants/constant.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class full_cell_solver{
public:
    Primitive primitive = Primitive::Line2;
    VectorXd point_coord;
    VectorXd result;
    const int iter = 10;
    const double tolerance = constant::tolerance;
    int an, ca;
    stiffness_separator sep;
    stiffness_anode anode;
    stiffness_cathode cathode;
    int step = 0;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    full_cell_solver(int an, int ca, const VectorXd& coord): point_coord(coord), an(an), ca(ca) {
        sep = stiffness_separator(coord, an, ca);
        anode = stiffness_anode(coord, an, ca);
        cathode = stiffness_cathode(coord, an, ca);
    }

    void calc(Eigen::Ref<MatrixXd>);
    void apply_boundary(Eigen::Ref<MatrixXd>, Eigen::SparseMatrix<double>&, Eigen::Ref<VectorXd>, bool);

};

#endif //FULL_CELL_SOLVER_H