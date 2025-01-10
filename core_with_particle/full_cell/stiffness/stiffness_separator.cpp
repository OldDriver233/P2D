#include "stiffness_separator.h"
#include "../../functions/functions.h"
#include "../../constants/constant.h"
#include <cmath>
#include <eigen3/Eigen/src/Core/Matrix.h>

template<bool use_temp>
void stiffness_separator::generate(const Eigen::Ref<MatrixXd> &u, 
                                   const Eigen::Ref<MatrixXd> &du, 
                                   const Eigen::Ref<MatrixXd> &c_s,
                                   std::vector<Eigen::Triplet<double>> &t, 
                                   Eigen::Ref<VectorXd> res, 
                                   bool is_first_step) {
    const int dim = 1, n = 2;
    long dof_cnt = this->surface_ca_coll - this->surface_an_coll + 1;
    long elem_cnt = dof_cnt - 1;
    int dof_cnt_eff = dof_cnt - (this->surface_ca_sep - this->surface_an_sep - 1);
    MatrixXd xs = get_integration_point<dim, n>();
    MatrixXd w = get_integration_weight<dim, n>();
    double dt = constant::dt;
    double eff_mat = std::pow(constant::epsilon_e_sep, constant::bruggeman);
    double d_ref = constant::de_sep;
    double d_eff = constant::de_sep / d_ref * eff_mat, epsilon = constant::epsilon_e_sep;
    double ce_int = constant::ce_int;
    double k_ref = constant::k_ref;
    double eff_1 = 1 / dt * constant::l_ref * constant::l_ref / d_ref;
    const double rho = 1100, cap = 700, lambda = 0.16;
    const int simd_size = surface_ca_sep - surface_an_sep + 1;
    std::vector<double> arr_kappa(2 * (simd_size - 1));
    std::vector<double> arr_d_kappa(2 * (simd_size - 1));

    if (!settings::use_customize_kappa) {
        for(int i = this->surface_an_sep - surface_an_coll; i < this->surface_ca_sep - surface_an_coll; ++i) {
            const int src_idx = i + surface_an_coll - surface_an_sep;
            MatrixXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0);
            //MatrixXd e_t = u({2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll, 2 * dof_cnt + 2 * dof_cnt_eff + i + 1 + this->surface_an_coll}, 0);
            for (int j = 0; j < n; j++) {
                const MatrixXd &N = cached_matrix_N[(i + surface_an_coll) * n + j];

                MatrixXd t_mat = N.transpose() * e_c;
                double ele_c_e = t_mat.sum();
                //t_mat = N.transpose() * e_t;
                //double ele_c_t = t_mat.sum();

                double k_eff = kappa(ele_c_e * ce_int) / constant::k_ref * eff_mat;
                double d_k_eff = d_kappa(ele_c_e * ce_int) * ce_int / constant::k_ref * eff_mat;
                arr_kappa[src_idx * n + j] = k_eff;
                arr_d_kappa[src_idx * n + j] = d_k_eff;
            }
        }
    } else {
        VectorXd kappa;
        VectorXd d_kappa;
        if constexpr(use_temp) {
            MatrixXd vars(2 * (simd_size - 1), 2);
            for(int i = this->surface_an_sep - surface_an_coll; i < this->surface_ca_sep - surface_an_coll; ++i) {
                const int src_idx = i + surface_an_coll - surface_an_sep;
                MatrixXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0);
                MatrixXd e_t = u({2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll, 2 * dof_cnt + 2 * dof_cnt_eff + i + 1 + this->surface_an_coll}, 0);
                for (int j = 0; j < n; j++) {
                    const MatrixXd &N = cached_matrix_N[(i + surface_an_coll) * n + j];

                    MatrixXd t_mat = N.transpose() * e_c;
                    double ele_c_e = t_mat.sum();
                    t_mat = N.transpose() * e_t;
                    double ele_c_t = t_mat.sum();

                    vars(src_idx * n + j, 0) = ele_c_e * ce_int;
                    vars(src_idx * n + j, 1) = ele_c_t;
                }
            }
            kappa = this->pfm->kappa.initial_node->eval(vars);
            d_kappa = this->pfm->kappa.initial_node->eval_deriv(vars, 0);
        } else {
            MatrixXd vars(2 * (simd_size - 1), 2);
            for (int i = 0; i < this->surface_an_sep - this->surface_an_coll; ++i) {
                VectorXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0);
                for (int j = 0; j < n; j++) {
                    const MatrixXd &N = cached_matrix_N[(i + surface_an_coll) * n + j];

                    MatrixXd t_mat = N.transpose() * e_c;
                    double ele_c_e = t_mat.sum();
                    vars(i * n + j, 0) = ele_c_e * ce_int;
                    vars(i * n + j, 1) = constant::T;
                }
            }
            kappa = this->pfm->kappa.initial_node->eval(vars);
            d_kappa = this->pfm->kappa.initial_node->eval_deriv(vars, 0);
        }
        for (int i = 0; i < 2 * (simd_size - 1); i++) {
            arr_kappa[i] = kappa(i) / constant::k_ref * eff_mat;
            arr_d_kappa[i] = d_kappa(i) * ce_int / constant::k_ref * eff_mat;
        }
    }

    for(int i = this->surface_an_sep - surface_an_coll; i < this->surface_ca_sep - surface_an_coll; ++i) {
        MatrixXd e_p = u({i, i + 1}, 0);
        MatrixXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0);
        MatrixXd e_t = MatrixXd::Ones(n, 1) * constant::T;
        MatrixXd e_dt = MatrixXd::Zero(n, 1);
        if constexpr(use_temp) {
            e_t = u({2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll, 2 * dof_cnt + 2 * dof_cnt_eff + i + 1 + this->surface_an_coll}, 0);
            e_dt = du({2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll, 2 * dof_cnt + 2 * dof_cnt_eff + i + 1 + this->surface_an_coll}, 0);
        }
        MatrixXd e_dc = du({dof_cnt + i, dof_cnt + i + 1}, 0);
        MatrixXd e_kpp = MatrixXd::Zero(n, n);
        MatrixXd e_kpc = MatrixXd::Zero(n, n);
        MatrixXd e_kcp = MatrixXd::Zero(n, n);
        MatrixXd e_kcc = MatrixXd::Zero(n, n);
        MatrixXd e_ktt = MatrixXd::Zero(n, n);
        MatrixXd e_ktp = MatrixXd::Zero(n, n);
        MatrixXd e_ktc = MatrixXd::Zero(n, n);
        MatrixXd e_rp = MatrixXd::Zero(n, 1);
        MatrixXd e_rc = MatrixXd::Zero(n, 1);
        MatrixXd e_rt = MatrixXd::Zero(n, 1);
        
        const int src_idx = i + surface_an_coll - surface_an_sep;

        for(int j = 0; j < n; ++j) {
            const MatrixXd &N = cached_matrix_N[(i + surface_an_coll) * n + j];
            const MatrixXd &dN = cached_matrix_dN[(i + surface_an_coll) * n + j];
            const MatrixXd &NdNT = cached_matrix_NdNT[(i + surface_an_coll) * n + j];
            const MatrixXd &NNT = cached_matrix_NNT[(i + surface_an_coll) * n + j];
            const MatrixXd &dNdNT = cached_matrix_dNdNT[(i + surface_an_coll) * n + j];
            const double det = cached_det_J[(i + surface_an_coll) * n + j];

            MatrixXd N_T = N.transpose();
            MatrixXd dN_T = dN.transpose();
            MatrixXd t_mat = N_T * e_c;
            double ele_c_e = t_mat.sum();
            t_mat = N.transpose() * e_t;
            double ele_c_t = t_mat.sum();

            double k_ref = constant::k_ref;
            double k_eff = arr_kappa[src_idx * n + j];
            double dk_dt = 2 * k_eff * constant::R / constant::F * (1 - constant::trans);
            double kd_eff = dk_dt * ele_c_t;
            double d_k_eff = arr_d_kappa[src_idx * n + j];
            double d_kd_eff = 2 * d_k_eff * constant::R * ele_c_t / constant::F * (1 - constant::trans);

            // phi part
            e_rp += k_eff * dNdNT * e_p * w(j) * det - kd_eff / ele_c_e * dNdNT * e_c * w(j) * det;
            e_kpp += k_eff * dNdNT * w(j) * det;
            e_kpc += d_k_eff * dNdNT * e_p * N_T * w(j) * det 
                     - kd_eff / ele_c_e * dNdNT * w(j) * det
                     - d_kd_eff / ele_c_e * dNdNT * e_c * N_T * w(j) * det
                     + kd_eff / (ele_c_e * ele_c_e) * dNdNT * e_c * N_T * w(j) * det;
            //e_kpp = MatrixXd::Identity(n, n);

            // c part
            e_rc += epsilon * eff_1 * NNT * e_dc * w(j) * det + d_eff * dNdNT * e_c * w(j) * det;
            e_kcc += epsilon * eff_1 * NNT * w(j) * det + d_eff * dNdNT * w(j) * det;
            //e_kcc = MatrixXd::Identity(n, n);

            // t part
            if constexpr(use_temp) {
                MatrixXd e_dp2 = (dN_T * e_p).array() * (dN_T * e_p).array();
                MatrixXd e_tdpdc = (N_T * e_t).array() * (dN_T * e_p).array() * (dN_T * e_c).array();
                MatrixXd e_dpdc = (dN_T * e_p).array() * (dN_T * e_c).array();
                MatrixXd e_tdc = (N_T * e_t).array() * (dN_T * e_c).array();
                MatrixXd e_tdp = (N_T * e_t).array() * (dN_T * e_p).array();
                e_ktt += rho * cap * constant::l_ref * constant::l_ref * NNT / constant::dt * w(j) * det
                         + lambda * dNdNT * w(j) * det;
                e_ktt += -(dk_dt * N * e_dpdc * N_T / ele_c_e);
                e_ktp += -(k_eff * N * 2 * dN_T * e_p * dN_T + dk_dt * N * e_tdc * dN_T / ele_c_e);
                e_ktc += -(dk_dt * N * e_tdp * dN_T / ele_c_e);
                e_rt += rho * cap * constant::l_ref * constant::l_ref * NNT * e_dt / constant::dt * w(j) * det
                        + lambda * dNdNT * e_t * w(j) * det;
                e_rt += -(k_eff * N * e_dp2 + dk_dt * N * e_tdpdc / ele_c_e) * k_ref * w(j) * det;
            }
        }

        {
            for(int j = 0; j < n; j++) {
                for(int l = 0; l < n; l++) {
                    //k(i + j, i + l) += e_kpp(j, l);
                    //k(i + j, i + l + dof_cnt) += e_kpc(j, l);
                    //k(i + j + dof_cnt, i + l) += e_kcp(j, l);
                    //k(i + j + dof_cnt, i + l + dof_cnt) += e_kcc(j, l);
                    t.emplace_back(i + j, i + l, e_kpp(j, l));
                    t.emplace_back(i + j, i + l + dof_cnt, e_kpc(j, l));
                    t.emplace_back(i + j + dof_cnt, i + l, e_kcp(j, l));
                    t.emplace_back(i + j + dof_cnt, i + l + dof_cnt, e_kcc(j, l));
                    if constexpr(use_temp) {
                        t.emplace_back(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j, 2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + l, e_ktt(j, l));
                        t.emplace_back(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j, i + l, e_ktp(j, l));
                        t.emplace_back(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j, i + l + dof_cnt, e_ktc(j, l));
                    }
                }
                res(i + j) += e_rp(j);
                res(i + j + dof_cnt) += e_rc(j);
                if constexpr(use_temp) res(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j) += e_rt(j);
            }
        }
    }
}

template void stiffness_separator::generate<true>(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                  bool);
template void stiffness_separator::generate<false>(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                  bool);