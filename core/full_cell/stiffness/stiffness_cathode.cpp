#include "stiffness_cathode.h"
#include "../../constants/constant.h"
#include "../../functions/functions.h"
#include <iostream>

void stiffness_cathode::generate(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> du, Eigen::Ref<MatrixXd> k, Eigen::Ref<VectorXd> res) {
    int dim = 1, n = 2;
    int dof_cnt = this->points.size();
    int dof_cnt_eff = dof_cnt - (this->surface_ca_sep - this->surface_an_sep - 1);
    int elem_cnt = this->points.size() - 1;
    MatrixXd xs = get_integration_point(dim, n);
    MatrixXd w = get_integration_weight(dim, n);
    double dt = constant::dt;
    double R_p = 2e-6, d_eff = 1.6478e-11, ds_eff = 1e-14, sigma_eff = 12.1173, epsilon = 0.385, epsilon_s = 0.59, k_eff = 0.02624, kd_eff = 8.5755e-4, a = 885000;
    const double F = 96485.3329;

    for(int i = this->surface_ca_sep; i < elem_cnt; i++) {
        int idx = i - this->surface_ca_sep + this->surface_an_sep + 1;
        MatrixXd e_p = u({i, i + 1}, 0); // phi_e
        MatrixXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0); // c_e
        MatrixXd e_s = u({2 * dof_cnt + idx, 2 * dof_cnt + idx + 1}, 0); // phi_s
        MatrixXd e_q = u({2 * dof_cnt + dof_cnt_eff + idx, 2 * dof_cnt + dof_cnt_eff + idx + 1}, 0); // j
        MatrixXd e_v = u({2 * dof_cnt + 2 * dof_cnt_eff + idx, 2 * dof_cnt + 2 * dof_cnt_eff + idx + 1}, 0); // c^avg
        MatrixXd e_a = u({2 * dof_cnt + 3 * dof_cnt_eff + idx, 2 * dof_cnt + 3 * dof_cnt_eff + idx + 1}, 0); // c^*

        MatrixXd e_dc = du({dof_cnt + i, dof_cnt + i + 1}, 0);
        MatrixXd e_dv = du({2 * dof_cnt + 2 * dof_cnt_eff + idx, 2 * dof_cnt + 2 * dof_cnt_eff + idx + 1}, 0);

        MatrixXd e_kvq = MatrixXd::Identity(n, n) * 3 / R_p;
        MatrixXd e_kvv = MatrixXd::Identity(n, n) / dt;
        MatrixXd e_kaa = MatrixXd::Identity(n, n);
        MatrixXd e_kav = MatrixXd::Identity(n, n) * -1;
        MatrixXd e_kaq = MatrixXd::Identity(n, n) * 0.2 * R_p / ds_eff;
        MatrixXd e_kqq = MatrixXd::Identity(n, n);

        MatrixXd e_kss = MatrixXd::Zero(n, n);
        MatrixXd e_ksq = MatrixXd::Zero(n, n);
        MatrixXd e_kcc = MatrixXd::Zero(n, n);
        MatrixXd e_kcq = MatrixXd::Zero(n, n);
        MatrixXd e_kpp = MatrixXd::Zero(n, n);
        MatrixXd e_kpc = MatrixXd::Zero(n, n);
        MatrixXd e_kpq = MatrixXd::Zero(n, n);
        MatrixXd e_kqp = MatrixXd::Zero(n, n);
        MatrixXd e_kqc = MatrixXd::Zero(n, n);
        MatrixXd e_kqs = MatrixXd::Zero(n, n);
        MatrixXd e_kqa = MatrixXd::Zero(n, n);

        MatrixXd e_rv = e_dv / dt + 3 * e_q / R_p;
        MatrixXd e_ra = e_a - e_v + e_q * 0.2 * R_p / d_eff;
        MatrixXd e_rs = MatrixXd::Zero(n, 1);
        MatrixXd e_rc = MatrixXd::Zero(n, 1);
        MatrixXd e_rq = MatrixXd::Zero(n, 1);
        MatrixXd e_rp = MatrixXd::Zero(n, 1);

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

            // s part
            e_kss += sigma_eff * dN * dN_T * w(j) * det;
            e_ksq += a * F * N * N_T * w(j) * det;
            e_rs += sigma_eff * dN * dN_T * e_s * w(j) * det + a * F * N * N_T * e_q * w(j) * det;

            // c part
            e_kcc += epsilon * N * N_T * w(j) * det / dt + d_eff * dN * dN_T * w(j) * det;
            e_kcq += -(1 - 0.4) * a * N * N_T * w(j) * det;
            e_rc += epsilon * N * N_T * e_dc * w(j) * det / dt + d_eff * dN * dN_T * e_c * w(j) * det - (1 - 0.4) * a * N * N_T * e_q * w(j) * det;

            // p part
            e_kpp += k_eff * dN * dN_T * w(j) * det;
            e_kpc += -kd_eff / lower * dN * dN_T * w(j) * det + kd_eff / (lower * lower) * dN * dN_T * e_c * N_T * w(j) * det;
            e_kpq += -a * F * N * N_T * w(j) * det;
            e_rp += k_eff * dN * dN_T * e_p * w(j) * det - kd_eff / lower * dN * dN_T * e_c * w(j) * det - a * F * N * N_T * e_q * w(j) * det;

            // q part
            double j0_v = j0(e_c(j), e_a(j), 2);
            double d_j0_a_v = d_j0_a(e_c(j), e_a(j), 2);
            double d_j0_e_v = d_j0_e(e_c(j), e_a(j), 2);
            double uoc_v = uoc(e_a(j) / 51554.0, 2);
            double d_uoc_v = d_uoc(e_a(j) / 51554.0, 2) / 51554.0;
            double bv_v = bv(e_s(j) - e_p(j) - uoc_v);
            double d_bv_v = d_bv(e_s(j) - e_p(j) - uoc_v);

            e_kqp(j, j) += -j0_v * d_bv_v * (-1);
            e_kqc(j, j) += -d_j0_e_v * bv_v;
            e_kqs(j, j) += -j0_v * d_bv_v;
            e_kqa(j, j) += -d_j0_a_v * bv_v - j0_v * d_bv_v * (-1) * d_uoc_v;
            e_rq(j) += e_q(j) - j0_v * bv_v;
        }

        for(int j = 0; j < n; j++) {
            for(int l = 0; l < n; l++) {
                k(i + j, i + l) += e_kpp(j, l);
                k(i + j, i + l + dof_cnt) += e_kpc(j, l);
                k(i + j, idx + l + 2 * dof_cnt + dof_cnt_eff) += e_kpq(j, l);

                k(i + j + dof_cnt, i + l + dof_cnt) += e_kcc(j, l);
                k(i + j + dof_cnt, idx + l + 2 * dof_cnt + dof_cnt_eff) += e_kcq(j, l); 

                k(idx + j + 2 * dof_cnt, idx + l + 2 * dof_cnt) += e_kss(j, l);
                k(idx + j + 2 * dof_cnt, idx + l + 2 * dof_cnt + dof_cnt_eff) += e_ksq(j, l);

                k(idx + j + 2 * dof_cnt + dof_cnt_eff, i + l) += e_kqp(j, l);
                k(idx + j + 2 * dof_cnt + dof_cnt_eff, i + l + dof_cnt) += e_kqc(j, l);
                k(idx + j + 2 * dof_cnt + dof_cnt_eff, idx + l + 2 * dof_cnt) += e_kqs(j, l);
                k(idx + j + 2 * dof_cnt + dof_cnt_eff, idx + l + 2 * dof_cnt + dof_cnt_eff) += e_kqq(j, l);
                k(idx + j + 2 * dof_cnt + dof_cnt_eff, idx + l + 2 * dof_cnt + 3 * dof_cnt_eff) += e_kqa(j, l);

                k(idx + j + 2 * dof_cnt + 2 * dof_cnt_eff, idx + l + 2 * dof_cnt + dof_cnt_eff) += e_kvq(j, l);
                k(idx + j + 2 * dof_cnt + 2 * dof_cnt_eff, idx + l + 2 * dof_cnt + 2 * dof_cnt_eff) += e_kvv(j, l);

                k(idx + j + 2 * dof_cnt + 3 * dof_cnt_eff, idx + l + 2 * dof_cnt + dof_cnt_eff) += e_kaq(j, l);
                k(idx + j + 2 * dof_cnt + 3 * dof_cnt_eff, idx + l + 2 * dof_cnt + 2 * dof_cnt_eff) += e_kav(j, l);
                k(idx + j + 2 * dof_cnt + 3 * dof_cnt_eff, idx + l + 2 * dof_cnt + 3 * dof_cnt_eff) += e_kaa(j, l);

                res(i + j) += e_rp(j);
                res(i + j + dof_cnt) += e_rc(j);
                res(idx + j + 2 * dof_cnt) += e_rs(j);
                res(idx + j + 2 * dof_cnt + dof_cnt_eff) += e_rq(j);
                res(idx + j + 2 * dof_cnt + 2 * dof_cnt_eff) += e_rv(j);
                res(idx + j + 2 * dof_cnt + 3 * dof_cnt_eff) += e_ra(j);
            }
        }

    }
}