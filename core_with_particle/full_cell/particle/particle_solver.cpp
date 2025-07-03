#include "particle_solver.h"

#include "../../io/settings/settings.h"
#include "../../step_control/StepControl.h"

void particle_solver::pre_calc(const Eigen::Ref<MatrixXd> &c_s) {
    this->last_cs = c_s;
}

template <bool constant_matrix>
double particle_solver::calc(Eigen::Ref<MatrixXd> c_s, const Eigen::Ref<MatrixXd> &u, int type, const mesh_reader& mesh, const dof_assigner& dof) {
    int particle_dof_size = constant::particle_segment + 1;
    double ret = 0.0;
    if (type == 1) {
        int i = 0;
        for (auto x: dof.particle_to_node) {
            if (mesh.anode_nodes.contains(x)) {
                if constexpr (!constant_matrix) {
                    double T;
                    if (settings::calc_temperature) {
                        T = u(dof.get_dof(x, 4), 0);
                    } else {
                        T = constant::t_ref;
                    }
                    double coeff = exp(constant::diffuse_energy_an / constant::R * (-1 / T + 1 / constant::t_ref));
                    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_A;
                    solver_A.compute(assembled_A_1 * coeff + assembled_A_2 / this->step_control->dt_now);
                    j_coeff = solver_A.solve(pre_j_coeff);
                    MatrixXd Bc_s = assembled_B / this->step_control->dt_now * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
                    c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) =
                    solver_A.solve(Bc_s) - j_coeff * u(dof.get_dof(x, 3), 0);
                    ret = j_coeff(constant::particle_segment);
                } else {
                    MatrixXd Bc_s = assembled_B / this->step_control->dt_now * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
                    c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) =
                    constant_solver.solve(Bc_s) - j_coeff * u(dof.get_dof(x, 3), 0);
                }
            }
            i++;
        }
    } else {
        int i = 0;
        for (auto x: dof.particle_to_node) {
            if (mesh.cathode_nodes.contains(x)) {
                if constexpr (!constant_matrix) {
                    double T;
                    if (settings::calc_temperature) {
                        T = u(dof.get_dof(x, 4), 0);
                    } else {
                        T = constant::t_ref;
                    }
                    double coeff = exp(constant::diffuse_energy_ca / constant::R * (-1 / T + 1 / constant::t_ref));
                    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_A;
                    solver_A.compute(assembled_A_1 * coeff + assembled_A_2 / this->step_control->dt_now);
                    j_coeff = solver_A.solve(pre_j_coeff);
                    MatrixXd Bc_s = assembled_B / this->step_control->dt_now * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
                    c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) =
                    solver_A.solve(Bc_s) - j_coeff * u(dof.get_dof(x, 3), 0);
                    ret = j_coeff(constant::particle_segment);
                } else {
                    MatrixXd Bc_s = assembled_B / this->step_control->dt_now * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
                    c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) =
                    constant_solver.solve(Bc_s) - j_coeff * u(dof.get_dof(x, 3), 0);
                }
            }
            i++;
        }
    }
    return ret;
}

template double particle_solver::calc<true>(Eigen::Ref<MatrixXd>, const Eigen::Ref<MatrixXd>&, int, const mesh_reader&, const dof_assigner&);
template double particle_solver::calc<false>(Eigen::Ref<MatrixXd>, const Eigen::Ref<MatrixXd>&, int, const mesh_reader&, const dof_assigner&);