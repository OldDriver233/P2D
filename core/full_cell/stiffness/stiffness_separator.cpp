#include "stiffness_separator.h"
#include "../../functions/functions.h"
#include "../../constants/constant.h"
#include <cmath>
#include <iostream>

void stiffness_separator::generate(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> du, Eigen::Ref<MatrixXd>k, Eigen::Ref<VectorXd> res) {
    int dim = 1, n = 2;
    long dof_cnt = this->points.size();
    long elem_cnt = this->points.size() - 1;
    MatrixXd xs = get_integration_point(dim, n);
    MatrixXd w = get_integration_weight(dim, n);
    double dt = constant::dt;
    double eff_mat = std::pow(constant::epsilon_e_sep, constant::bruggeman);
    double d_ref = constant::de_sep;
    double d_eff = constant::de_sep / d_ref * eff_mat, epsilon = constant::epsilon_e_sep;
    double ce_int = constant::ce_int;

    for(int i = this->surface_an_sep; i < this->surface_ca_sep; i++) {
        MatrixXd e_p = u({i, i + 1}, 0);
        MatrixXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0);
        MatrixXd e_dc = du({dof_cnt + i, dof_cnt + i + 1}, 0);
        MatrixXd e_kpp = MatrixXd::Zero(n, n);
        MatrixXd e_kpc = MatrixXd::Zero(n, n);
        MatrixXd e_kcp = MatrixXd::Zero(n, n);
        MatrixXd e_kcc = MatrixXd::Zero(n, n);
        MatrixXd e_rp = MatrixXd::Zero(n, 1);
        MatrixXd e_rc = MatrixXd::Zero(n, 1);

        for(int j = 0; j < n; j++) {
            const MatrixXd &N = cached_matrix_N[i * n + j];
            const MatrixXd &dN = cached_matrix_dN[i * n + j];
            MatrixXd N_T = N.transpose();
            MatrixXd dN_T = dN.transpose();
            double det = cached_det_J[i * n + j];
            double s = xs(j);
            MatrixXd t_mat = N_T * e_c;
            // t_mat should be 1x1
            double lower = t_mat.sum();
            double ce_int = constant::ce_int;

            double k_ref = constant::k_ref;
            double k_eff = kappa(lower * ce_int) / k_ref * eff_mat, kd_eff = 2 * k_eff * constant::R * constant::T / constant::F * (1 - constant::trans);
            double d_k_eff = d_kappa(lower * ce_int) * ce_int / k_ref * eff_mat;
            double d_kd_eff = 2 * d_k_eff * constant::R * constant::T / constant::F * (1 - constant::trans);
            
            // phi part
            e_rp += k_eff * dN * dN_T * e_p * w(j) * det - kd_eff / lower * dN * dN_T * e_c * w(j) * det;
            e_kpp += k_eff * dN * dN_T * w(j) * det;
            e_kpc += d_k_eff * dN * dN_T * e_p * N_T * w(j) * det 
                     - kd_eff / lower * dN * dN_T * w(j) * det
                     - d_kd_eff / lower * dN * dN_T * e_c * N_T * w(j) * det
                     + kd_eff / (lower * lower) * dN * dN_T * e_c * N_T * w(j) * det;
            //e_kpp = MatrixXd::Identity(n, n);

            // c part
            double eff_1 = 1 / dt * constant::l_ref * constant::l_ref / d_ref;
            e_rc += epsilon * eff_1 * N * N_T * e_dc * w(j) * det + d_eff * dN * dN_T * e_c * w(j) * det;
            e_kcc += epsilon * eff_1 * N * N_T * w(j) * det + d_eff * dN * dN_T * w(j) * det;
            //e_kcc = MatrixXd::Identity(n, n);
        }

        for(int j = 0; j < n; j++) {
            for(int l = 0; l < n; l++) {
                k(i + j, i + l) += e_kpp(j, l);
                k(i + j, i + l + dof_cnt) += e_kpc(j, l);
                k(i + j + dof_cnt, i + l) += e_kcp(j, l);
                k(i + j + dof_cnt, i + l + dof_cnt) += e_kcc(j, l);
            }
            res(i + j) += e_rp(j);
            res(i + j + dof_cnt) += e_rc(j);
        }
    }
}