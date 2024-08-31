#include "stiffness_separator.h"
#include "../../constants/constant.h"
#include <iostream>

void stiffness_separator::generate(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> du, Eigen::Ref<MatrixXd>k, Eigen::Ref<VectorXd> res) {
    int dim = 1, n = 2;
    long dof_cnt = this->points.size();
    long elem_cnt = this->points.size() - 1;
    MatrixXd xs = get_integration_point(dim, n);
    MatrixXd w = get_integration_weight(dim, n);
    double dt = constant::dt;
    double k_eff = 0.2988, kd_eff = 9.053e-3, d_eff = 7.5e-10, epsilon = 0.724;

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
            
            // phi part
            e_rp += k_eff * dN * dN_T * e_p * w(j) * det - kd_eff / lower * dN * dN_T * e_c * w(j) * det;
            e_kpp += k_eff * dN * dN_T * w(j) * det;
            e_kpc += -kd_eff / lower * dN * dN_T * w(j) * det + kd_eff / (lower * lower) * dN * dN_T * e_c * N_T * w(j) * det;

            // c part
            e_rc += epsilon * N * N_T * e_dc * w(j) * det / dt + d_eff * dN * dN_T * e_c * w(j) * det;
            e_kcc += epsilon * N * N_T * w(j) * det / dt + d_eff * dN * dN_T * w(j) * det;
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