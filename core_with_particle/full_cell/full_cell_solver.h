#ifndef FULL_CELL_SOLVER_H
#define FULL_CELL_SOLVER_H
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include "../constants/constant.h"
#include "../shaping/primitive_type.h"
#include "../mesh/mesh_reader.h"
#include "../mesh/dof_assigner.h"
#include "stiffness/pre_calc.h"
#include "stiffness/stiffness_anode.h"
#include "stiffness/stiffness_cathode.h"
#include "stiffness/stiffness_separator.h"
#include "stiffness/stiffness_anode_collector.h"
#include "stiffness/stiffness_cathode_collector.h"
#include "particle/particle_solver.h"
#include "../functions/function_manager.h"
#include "../io/output/output_manager.h"
#include "../step_control/StepControl.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class full_cell_solver {
public:
    const mesh_reader& mesh;
    const dof_assigner& dof;
    pre_calc_shapes shapes;
    const int iter = 10;
    const double tolerance = constant::tolerance;
    //int an, ca, ancoll, cacoll;
    particle_solver anode_particle, cathode_particle;
    stiffness_separator sep;
    stiffness_anode anode;
    stiffness_cathode cathode;
    stiffness_anode_collector anode_collector;
    stiffness_cathode_collector cathode_collector;
    int step = 0;
    double current_time = 0.0;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    FunctionManager manager;
    //output_manager Q_ohm, Q_rxn, Q_rev;
    //output_manager eta;
    std::vector<double> temp;
    StepControl* step_control;


    full_cell_solver(const mesh_reader& mesh, const dof_assigner& dof, VectorXd &particle_coord, StepControl* st)
        : mesh(mesh), dof(dof), shapes(mesh, dof),
        anode_particle(particle_coord, constant::ds_an, constant::c_max_an, st),
        cathode_particle(particle_coord, constant::ds_ca, constant::c_max_ca, st),
        anode(mesh, dof, shapes, -anode_particle.j_coeff(constant::particle_segment), &manager, st),
        cathode(mesh, dof, shapes, -cathode_particle.j_coeff(constant::particle_segment), &manager, st),
        sep(mesh, dof, shapes, &manager, st),
        anode_collector(mesh, dof, shapes, st),
        cathode_collector(mesh, dof, shapes, st),
        step_control(st)
    {
        /*
        element_coord = VectorXd::Zero(cacoll - ancoll);
        for (auto i = 0; i < element_coord.size(); i++) {
            element_coord(i) = (point_coord(i + ancoll) + point_coord(i + ancoll + 1)) / 2 - point_coord(ancoll);
        }
        sep = stiffness_separator(coord, an, ca, ancoll, cacoll, &manager, st);
        anode = stiffness_anode(coord, an, ca, ancoll, cacoll, -anode_particle.j_coeff(constant::particle_segment), &manager, st);
        cathode = stiffness_cathode(coord, an, ca, ancoll, cacoll, -cathode_particle.j_coeff(constant::particle_segment), &manager, st);
        anode_collector = stiffness_anode_collector(coord, an, ca, ancoll, cacoll, st);
        cathode_collector = stiffness_cathode_collector(coord, an, ca, ancoll, cacoll, st);
        Q_ohm = output_manager(element_coord);
        Q_rxn = output_manager(element_coord);
        Q_rev = output_manager(element_coord);
        eta = output_manager(element_coord);
        */
    }

    void print_detail();
    void calc(Eigen::Ref<MatrixXd>, Eigen::Ref<MatrixXd>, double, bool);
    void apply_boundary(Eigen::Ref<MatrixXd>, Eigen::SparseMatrix<double> &,
                        Eigen::Ref<VectorXd>, bool);
};

#endif // FULL_CELL_SOLVER_H
