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
                                   bool is_first_step,
                                   std::vector<double>& temp) {
    const int dim = 1, n = 2;
    long dof_cnt = this->surface_ca_coll - this->surface_an_coll + 1;
    long elem_cnt = dof_cnt - 1;
    int dof_cnt_eff = dof_cnt - (this->surface_ca_sep - this->surface_an_sep - 1);
    MatrixXd xs = get_integration_point<dim, n>();
    MatrixXd w = get_integration_weight<dim, n>();
    double dt = constant::dt;
    double eff_mat = std::pow(constant::epsilon_e_sep, constant::bruggeman);
    double d_ref = constant::de_sep;
    //double d_eff = constant::de_sep / d_ref * eff_mat;
    double epsilon = constant::epsilon_e_sep;
    double ce_int = constant::ce_int;
    double k_ref = constant::k_ref;
    double eff_1 = 1 / dt * constant::l_ref * constant::l_ref / d_ref;
    const double rho = constant::density_sep, cap = constant::capacity_sep, lambda = constant::lambda_sep;
    const int simd_size = surface_ca_sep - surface_an_sep + 1;
    std::vector<double> arr_kappa(2 * (simd_size - 1));
    std::vector<double> arr_d_kappa(2 * (simd_size - 1));
    std::vector<double> arr_d_eff(2 * (simd_size - 1));

    MatrixXd vars(2 * (simd_size - 1), 2);
    if constexpr(use_temp) {
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
    } else {
        for(int i = this->surface_an_sep - surface_an_coll; i < this->surface_ca_sep - surface_an_coll; ++i) {
            const int src_idx = i + surface_an_coll - surface_an_sep;
            VectorXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0);
            for (int j = 0; j < n; j++) {
                const MatrixXd &N = cached_matrix_N[(i + surface_an_coll) * n + j];

                MatrixXd t_mat = N.transpose() * e_c;
                double ele_c_e = t_mat.sum();
                vars(src_idx * n + j, 0) = ele_c_e * ce_int;
                vars(src_idx * n + j, 1) = constant::T;
            }
        }
    }

    if (!settings::use_customize_kappa) {
        for(int i = this->surface_an_sep - surface_an_coll; i < this->surface_ca_sep - surface_an_coll; ++i) {
            const int src_idx = i + surface_an_coll - surface_an_sep;
            for (int j = 0; j < n; j++) {
                double k_eff = kappa(vars(src_idx * n + j)) / constant::k_ref * eff_mat;
                double d_k_eff = d_kappa(vars(src_idx * n + j)) * ce_int / constant::k_ref * eff_mat;
                arr_kappa[src_idx * n + j] = k_eff;
                arr_d_kappa[src_idx * n + j] = d_k_eff;
            }
        }
    } else {
        VectorXd kappa = this->pfm->kappa.initial_node->eval(vars);
        VectorXd d_kappa = this->pfm->kappa.initial_node->eval_deriv(vars, 0);
        for (int i = 0; i < 2 * (simd_size - 1); i++) {
            arr_kappa[i] = kappa(i) / constant::k_ref * eff_mat;
            arr_d_kappa[i] = d_kappa(i) * ce_int / constant::k_ref * eff_mat;
        }
    }
    if (!settings::use_customize_diffuse) {
        for (int i = 0; i < 2 * (simd_size - 1); i++) {
            arr_d_eff[i] = constant::de_sep / d_ref * eff_mat;
        }
    } else {
        VectorXd d_l = this->pfm->electrolyte_diffuse.initial_node->eval(vars);
        for (int i = 0; i < 2 * (simd_size - 1); i++) {
            arr_d_eff[i] = d_l(i) / d_ref * eff_mat;
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
            e_rc += epsilon * eff_1 * NNT * e_dc * w(j) * det + arr_d_eff[src_idx * n + j] * dNdNT * e_c * w(j) * det;
            e_kcc += epsilon * eff_1 * NNT * w(j) * det + arr_d_eff[src_idx * n + j] * dNdNT * w(j) * det;
            //e_kcc = MatrixXd::Identity(n, n);

            // t part
            if constexpr(use_temp) {
                double dNe_p = (dN_T * e_p).sum();
                double Ne_t = (N_T * e_t).sum();
                double dNe_c = (dN_T * e_c).sum();
                double e_dp2 = dNe_p * dNe_p;
                double e_tdpdc = Ne_t * dNe_p * dNe_c;
                double e_dpdc = dNe_p * dNe_c;
                double e_tdc = Ne_t * dNe_c;
                double e_tdp = Ne_t * dNe_p;
                e_ktt += rho * cap * constant::l_ref * constant::l_ref * NNT / constant::dt * w(j) * det
                         + lambda * dNdNT * w(j) * det;
                e_rt += rho * cap * constant::l_ref * constant::l_ref * NNT * e_dt / constant::dt * w(j) * det
                         + lambda * dNdNT * e_t * w(j) * det;
                if (!is_first_step) {
                    e_ktt += (dk_dt * N * e_dpdc * N_T / ele_c_e) * k_ref;
                    e_ktp += -(k_eff * N * 2 * dN_T * e_p * dN_T - dk_dt * N * e_tdc * dN_T / ele_c_e) * k_ref;
                    e_ktc += (dk_dt * N * e_tdp * dN_T / ele_c_e) * k_ref;
                    e_rt += -(k_eff * N * e_dp2 - dk_dt * N * e_tdpdc / ele_c_e) * k_ref * w(j) * det;
                }
                if (j == 0) {
                    temp.push_back(((k_eff * e_dp2 - dk_dt * e_tdpdc / ele_c_e) * k_ref) / (constant::l_ref * constant::l_ref));
                    temp.push_back(0);
                    temp.push_back(0);
                    temp.push_back(0);
                }
            } else if (j == 0) {
                temp.push_back(0);
                temp.push_back(0);
                temp.push_back(0);
                temp.push_back(0);
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
                    //t.emplace_back(i + j + dof_cnt, i + l, e_kcp(j, l));
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
                  bool, std::vector<double> &);
template void stiffness_separator::generate<false>(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                  bool, std::vector<double> &);