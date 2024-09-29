#include "stiffness_cathode.h"
#include "../../constants/constant.h"
#include "../../functions/functions.h"
#include <iostream>

void stiffness_cathode::generate(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> du, std::vector<Eigen::Triplet<double>> &t, Eigen::Ref<VectorXd> res, bool is_first_step) {
    int dim = 1, n = 2;
    int dof_cnt = this->points.size();
    int dof_cnt_eff = dof_cnt - (this->surface_ca_sep - this->surface_an_sep - 1);
    int elem_cnt = this->points.size() - 1;
    MatrixXd xs = get_integration_point(dim, n);
    MatrixXd w = get_integration_weight(dim, n);
    double dt = constant::dt;
    double R_p = constant::r_p;
    double a = 3 * constant::epsilon_s_ca / constant::r_p;
    double eff_mat = std::pow(constant::epsilon_e_ca, constant::bruggeman);
    double eff_mat_s = std::pow(constant::epsilon_s_ca, constant::bruggeman);
    double d_ref = constant::de_ca;
    double d_eff = constant::de_ca / d_ref * eff_mat;
    double ds_eff = constant::ds_ca / d_ref;
    double sigma_ref = constant::sigma_ca * eff_mat_s;
    double sigma_eff = constant::sigma_ca / sigma_ref * eff_mat_s;
    double epsilon = constant::epsilon_e_ca;
    double epsilon_s = constant::epsilon_s_ca;
    double c_max = constant::c_max_ca;
    double ce_int = constant::ce_int;
    double j_ref = constant::j_ref;
    const double F = constant::F;

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

        MatrixXd e_kvq = MatrixXd::Identity(n, n) * 3 / R_p / c_max * j_ref;
        MatrixXd e_kvv = MatrixXd::Identity(n, n) / dt;
        MatrixXd e_kaa = MatrixXd::Identity(n, n);
        MatrixXd e_kav = MatrixXd::Identity(n, n) * -1;
        MatrixXd e_kaq = MatrixXd::Identity(n, n) * 0.2 * R_p / ds_eff / c_max * j_ref / d_ref;
        //MatrixXd e_kaq = MatrixXd::Zero(n, n);
        MatrixXd e_kqq = MatrixXd::Identity(n, n);
        //MatrixXd e_kvq = MatrixXd::Zero(n, n);
        //MatrixXd e_kvv = MatrixXd::Identity(n, n);
        //MatrixXd e_kaa = MatrixXd::Identity(n, n);
        //MatrixXd e_kav = MatrixXd::Zero(n, n);
        //MatrixXd e_kaq = MatrixXd::Zero(n, n);
        //MatrixXd e_kqq = MatrixXd::Identity(n, n);

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

        MatrixXd e_rv = e_dv / dt + 3 * e_q / R_p / c_max * j_ref;
        MatrixXd e_ra = e_a - e_v + e_q * 0.2 * R_p / ds_eff / c_max * j_ref / d_ref;
        //MatrixXd e_ra = MatrixXd::Zero(n, 1);
        //MatrixXd e_rv = MatrixXd::Zero(n, 1);
        //MatrixXd e_ra = MatrixXd::Zero(n, 1);
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

            double k_ref = constant::k_ref;
            double k_eff = kappa(lower * ce_int) / k_ref * eff_mat, kd_eff = 2 * k_eff * constant::R * constant::T / constant::F * (1 - constant::trans);
            double d_k_eff = d_kappa(lower * ce_int) * ce_int / k_ref * eff_mat;
            double d_kd_eff = 2 * d_k_eff * constant::R * constant::T / constant::F * (1 - constant::trans);

            // s part
            double eff_1 = a * F * constant::l_ref * constant::l_ref / sigma_ref;
            e_kss += sigma_eff * dN * dN_T * w(j) * det;
            e_ksq += eff_1 * N * N_T * w(j) * det * j_ref;
            e_rs += sigma_eff * dN * dN_T * e_s * w(j) * det + eff_1 * N * N_T * e_q * w(j) * det * j_ref;
            //e_kss = MatrixXd::Identity(2, 2);

            // c part
            double eff_2 = a * constant::l_ref * constant::l_ref * (1 - constant::trans) / d_ref;
            double eff_3 = 1 / dt * constant::l_ref * constant::l_ref / d_ref;
            e_kcc += epsilon * eff_3 * N * N_T * w(j) * det + d_eff * dN * dN_T * w(j) * det;
            e_kcq += -eff_2 * N * N_T * w(j) * det * j_ref / ce_int;
            e_rc += epsilon * eff_3 * N * N_T * e_dc * w(j) * det + d_eff * dN * dN_T * e_c * w(j) * det - eff_2 * N * N_T * e_q * w(j) * det * j_ref / ce_int;
            //e_kcc = MatrixXd::Identity(2, 2);

            // p part
            double eff_4 = a * F * constant::l_ref * constant::l_ref / k_ref;
            e_kpp += k_eff * dN * dN_T * w(j) * det;
            e_kpc += d_k_eff * dN * dN_T * e_p * N_T * w(j) * det 
                     - kd_eff / lower * dN * dN_T * w(j) * det
                     - d_kd_eff / lower * dN * dN_T * e_c * N_T * w(j) * det
                     + kd_eff / (lower * lower) * dN * dN_T * e_c * N_T * w(j) * det;
            e_kpq += -eff_4 * N * N_T * w(j) * det * j_ref;
            e_rp += k_eff * dN * dN_T * e_p * w(j) * det - kd_eff / lower * dN * dN_T * e_c * w(j) * det - eff_4 * N * N_T * e_q * w(j) * det * j_ref;
            //e_kpp = MatrixXd::Identity(2, 2);

            // q part
            double j0_v = j0(e_c(j), e_a(j), 2);
            double d_j0_a_v = d_j0_a(e_c(j), e_a(j), 2);
            double d_j0_e_v = d_j0_e(e_c(j), e_a(j), 2);
            double uoc_v = uoc(e_a(j), 2);
            double d_uoc_v = d_uoc(e_a(j), 2);
            double bv_v = bv(e_s(j) - e_p(j) - uoc_v);
            double d_bv_v = d_bv(e_s(j) - e_p(j) - uoc_v);
            double ce_root = std::sqrt(ce_int);

            e_kqp(j, j) += -j0_v * d_bv_v * (-1) * c_max * ce_root / j_ref;
            e_kqc(j, j) += -d_j0_e_v * bv_v * c_max * ce_root / j_ref;
            e_kqs(j, j) += -j0_v * d_bv_v * c_max * ce_root / j_ref;
            e_kqa(j, j) += -d_j0_a_v * bv_v * c_max * ce_root / j_ref - j0_v * d_bv_v * (-1) * d_uoc_v * c_max * ce_root / j_ref;
            e_rq(j) += e_q(j) - j0_v * bv_v * c_max * ce_root / j_ref;
            //e_rq(j) += e_q(j);
            //std::cout<<e_c(j)<<" "<<e_a(j)<<" "<<j0_v<<" "<<d_j0_a_v<<" "<<d_j0_e_v<<" "<<uoc_v<<" "<<d_uoc_v<<" "<<bv_v<<" "<<d_bv_v<<" "<<std::endl;
        }

        for(int j = 0; j < n; j++) {
            for(int l = 0; l < n; l++) {
                /*
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
                */
                t.push_back(Eigen::Triplet<double>(i + j, i + l, e_kpp(j, l)));
                t.push_back(Eigen::Triplet<double>(i + j, i + l + dof_cnt, e_kpc(j, l)));
                t.push_back(Eigen::Triplet<double>(i + j, idx + l + 2 * dof_cnt + dof_cnt_eff, e_kpq(j, l)));

                t.push_back(Eigen::Triplet<double>(i + j + dof_cnt, i + l + dof_cnt, e_kcc(j, l)));
                t.push_back(Eigen::Triplet<double>(i + j + dof_cnt, idx + l + 2 * dof_cnt + dof_cnt_eff, e_kcq(j, l)));

                if((!is_first_step) || i != elem_cnt - 1 || j != n - 1) {
                    t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt, idx + l + 2 * dof_cnt, e_kss(j, l)));
                    t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt, idx + l + 2 * dof_cnt + dof_cnt_eff, e_ksq(j, l)));
                }

                t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt + dof_cnt_eff, i + l, e_kqp(j, l)));
                t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt + dof_cnt_eff, i + l + dof_cnt, e_kqc(j, l)));
                t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt + dof_cnt_eff, idx + l + 2 * dof_cnt, e_kqs(j, l)));
                t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt + dof_cnt_eff, idx + l + 2 * dof_cnt + dof_cnt_eff, e_kqq(j, l)));
                t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt + dof_cnt_eff, idx + l + 2 * dof_cnt + 3 * dof_cnt_eff, e_kqa(j, l)));

                t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt + 2 * dof_cnt_eff, idx + l + 2 * dof_cnt + dof_cnt_eff, e_kvq(j, l)));
                t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt + 2 * dof_cnt_eff, idx + l + 2 * dof_cnt + 2 * dof_cnt_eff, e_kvv(j, l)));

                t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt + 3 * dof_cnt_eff, idx + l + 2 * dof_cnt + dof_cnt_eff, e_kaq(j, l)));
                t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt + 3 * dof_cnt_eff, idx + l + 2 * dof_cnt + 2 * dof_cnt_eff, e_kav(j, l)));
                t.push_back(Eigen::Triplet<double>(idx + j + 2 * dof_cnt + 3 * dof_cnt_eff, idx + l + 2 * dof_cnt + 3 * dof_cnt_eff, e_kaa(j, l)));
            }
            res(i + j) += e_rp(j);
            res(i + j + dof_cnt) += e_rc(j);
            res(idx + j + 2 * dof_cnt) += e_rs(j);
            res(idx + j + 2 * dof_cnt + dof_cnt_eff) += e_rq(j);
            res(idx + j + 2 * dof_cnt + 2 * dof_cnt_eff) += e_rv(j);
            res(idx + j + 2 * dof_cnt + 3 * dof_cnt_eff) += e_ra(j);
        }

    }
}