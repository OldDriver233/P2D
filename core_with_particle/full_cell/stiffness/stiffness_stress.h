#ifndef STIFFNESS_STRESS_H
#define STIFFNESS_STRESS_H
#include "../../functions/function_manager.h"
#include "../../mesh/mesh_reader.h"
#include "../../mesh/dof_assigner.h"
#include "pre_calc.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_stress {
public:
    const mesh_reader &mesh;
    const dof_assigner &dof;
    const pre_calc_shapes &shapes;

    stiffness_stress(const mesh_reader &mesh, const dof_assigner &dof,
                     const pre_calc_shapes &shapes): mesh(mesh), dof(dof), shapes(shapes) {
    }

    void generate(std::vector<Eigen::Triplet<double> > &t, std::vector<Eigen::Triplet<double> > &l);

    void generate_residue(const Eigen::Ref<MatrixXd> &u, const std::vector<double> avg_c_an,
                          const std::vector<double> avg_c_ca, Eigen::Ref<VectorXd> res,
                          const Eigen::Ref<Eigen::SparseMatrix<double>>& K);
};

#endif //STIFFNESS_STRESS_H
