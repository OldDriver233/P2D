#include "particle_solver.h"

void particle_solver::pre_calc(const Eigen::Ref<MatrixXd> &c_s) {
    this->last_cs = c_s;
}

void particle_solver::calc(Eigen::Ref<MatrixXd> c_s, const Eigen::Ref<MatrixXd> &u, int pt_size, int an, int ca, int type) {
    int eff_size = pt_size - (ca - an - 1);
    int particle_dof_size = constant::particle_segment + 1;
    if(type == 1) {
        for(int i = 0; i <= an; i++) {
            MatrixXd Bc_s = assembled_B * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
            c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) = 
            solver_A.solve(Bc_s) - j_coeff * u(2 * pt_size + eff_size + i, 0);
        }
    } else {
        for(int i = ca; i < pt_size; i++) {
            int idx = i - ca + an + 1;
            MatrixXd Bc_s = assembled_B * last_cs.block(idx * particle_dof_size, 0, particle_dof_size, 1);
            c_s.block(idx * particle_dof_size, 0, particle_dof_size, 1) = 
            solver_A.solve(Bc_s) - j_coeff * u(2 * pt_size + eff_size + idx, 0);
        }
    }
}