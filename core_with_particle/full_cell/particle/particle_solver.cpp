#include "particle_solver.h"

void particle_solver::pre_calc(const Eigen::Ref<MatrixXd> &c_s) {
    this->last_cs = c_s;
}

template <bool constant_matrix>
void particle_solver::calc(Eigen::Ref<MatrixXd> c_s, const Eigen::Ref<MatrixXd> &u, int pt_size, int an, int ca, int type, int temp) {
    int eff_size = pt_size - (ca - an - 1);
    int particle_dof_size = constant::particle_segment + 1;
    if(type == 1) {
        for(int i = 0; i <= an; i++) {
            if constexpr (!constant_matrix) {
                double coeff = exp(5000 / constant::R * (1 / u(temp + i, 0) - 1 / constant::t_ref));
                Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_A;
                solver_A.compute(assembled_A_1 * coeff + assembled_A_2);
                j_coeff = solver_A.solve(pre_j_coeff);
                MatrixXd Bc_s = assembled_B * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
                c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) =
                solver_A.solve(Bc_s) - j_coeff * u(2 * pt_size + eff_size + i, 0);
            } else {
                MatrixXd Bc_s = assembled_B * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
                c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) =
                constant_solver.solve(Bc_s) - j_coeff * u(2 * pt_size + eff_size + i, 0);
            }
        }
    } else {
        for(int i = ca; i < pt_size; i++) {
            if constexpr (!constant_matrix) {
                double coeff = exp(5000 / constant::R * (1 / u(temp + i, 0) - 1 / constant::t_ref));
                Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_A;
                solver_A.compute(assembled_A_1 * coeff + assembled_A_2);
                j_coeff = solver_A.solve(pre_j_coeff);
                int idx = i - ca + an + 1;
                MatrixXd Bc_s = assembled_B * last_cs.block(idx * particle_dof_size, 0, particle_dof_size, 1);
                c_s.block(idx * particle_dof_size, 0, particle_dof_size, 1) =
                solver_A.solve(Bc_s) - j_coeff * u(2 * pt_size + eff_size + idx, 0);
            } else {
                int idx = i - ca + an + 1;
                MatrixXd Bc_s = assembled_B * last_cs.block(idx * particle_dof_size, 0, particle_dof_size, 1);
                c_s.block(idx * particle_dof_size, 0, particle_dof_size, 1) =
                constant_solver.solve(Bc_s) - j_coeff * u(2 * pt_size + eff_size + idx, 0);
            }
        }
    }
}

template void particle_solver::calc<true>(Eigen::Ref<MatrixXd>, const Eigen::Ref<MatrixXd>&, int, int, int, int, int);
template void particle_solver::calc<false>(Eigen::Ref<MatrixXd>, const Eigen::Ref<MatrixXd>&, int, int, int, int, int);