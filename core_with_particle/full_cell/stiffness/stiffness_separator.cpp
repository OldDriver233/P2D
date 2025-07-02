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
    const int dim = get_dim(mesh.p_type), n = get_nodes(mesh.p_type);
    std::size_t node_cnt = mesh.node_count;

    MatrixXd xs;
    MatrixXd w;

    if (dim == 1) {
        xs = get_integration_point<1, 2>();
        w = get_integration_weight<1, 2>();
    } else if (dim == 2) {
        if (n == 3) {
            xs = get_integration_point<2, 3>();
            w = get_integration_weight<2, 3>();
        } else {
            xs = get_integration_point<2, 4>();
            w = get_integration_weight<2, 4>();
        }
    }

    double dt = this->st->dt_now;
    double eff_mat = std::pow(constant::epsilon_e_sep, constant::bruggeman);
    double d_ref = constant::de_sep;
    //double d_eff = constant::de_sep / d_ref * eff_mat;
    double epsilon = constant::epsilon_e_sep;
    double ce_int = constant::ce_int;
    double k_ref = constant::k_ref;
    double eff_1 = 1 / dt * constant::l_ref * constant::l_ref / d_ref;
    const double rho = constant::density_sep, cap = constant::capacity_sep, lambda = constant::lambda_sep;

    std::vector<double> arr_kappa(n * mesh.separator_elements.size());
    std::vector<double> arr_d_kappa(n * mesh.separator_elements.size());
    std::vector<double> arr_d_eff(n * mesh.separator_elements.size());

    MatrixXd vars(n * mesh.separator_elements.size(), 2);
    if constexpr(use_temp) {
        int i = 0;
        for (auto e: mesh.separator_elements) {
            MatrixXd e_c(n, 1);
            MatrixXd e_t(n, 1);

            for (int j = 0; j < n; j++) {
                int node_id = mesh.elements[e * n + j];
                e_c(j) = u(dof.get_dof(node_id, 1), 0);
                e_t(j) = u(dof.get_dof(node_id, 4), 0);
            }
            for (int j = 0; j < n; j++) {
                const MatrixXd &N = shapes.cached_matrix_N[e * n + j];
                MatrixXd t_mat = N.transpose() * e_c;
                double ele_c_e = t_mat.sum();
                t_mat = N.transpose() * e_t;
                double ele_c_t = t_mat.sum();
                vars(i * n + j, 0) = ele_c_e * ce_int;
                vars(i * n + j, 1) = ele_c_t;
            }
            i++;
        }
    } else {
        int i = 0;
        for (auto e: mesh.separator_elements) {
            MatrixXd e_c(n, 1);
            for (int j = 0; j < n; j++) {
                int node_id = mesh.elements[e * n + j];
                e_c(j) = u(dof.get_dof(node_id, 1), 0);
            }
            for (int j = 0; j < n; j++) {
                const MatrixXd &N = shapes.cached_matrix_N[e * n + j];
                MatrixXd t_mat = N.transpose() * e_c;
                double ele_c_e = t_mat.sum();
                vars(i * n + j, 0) = ele_c_e * ce_int;
                vars(i * n + j, 1) = constant::T;
            }
            i++;
        }
    }

    if (!settings::use_customize_kappa) {
        for (int i = 0; i < mesh.separator_elements.size(); i++) {
            for (int j = 0; j < n; j++) {
                double k_eff = kappa(vars(i * n + j)) / constant::k_ref * eff_mat;
                double d_k_eff = d_kappa(vars(i * n + j)) * ce_int / constant::k_ref * eff_mat;
                arr_kappa[i * n + j] = k_eff;
                arr_d_kappa[i * n + j] = d_k_eff;
            }
        }
    } else {
        for (int i = 0; i < mesh.separator_elements.size(); i++) {
            for (int j = 0; j < n; j++) {
                double k_eff = pfm->f_kappa(vars(i * n + j, 0), vars(i * n + j, 1)) / constant::k_ref * eff_mat;
                double d_k_eff = pfm->f_d_kappa(vars(i * n + j, 0), vars(i * n + j, 1)) * ce_int / constant::k_ref * eff_mat;
                arr_kappa[i * n + j] = k_eff;
                arr_d_kappa[i * n + j] = d_k_eff;
            }
        }
    }

    if (!settings::use_customize_diffuse) {
        for (int i = 0; i < n * mesh.separator_elements.size(); i++) {
            arr_d_eff[i] = constant::de_sep / d_ref * eff_mat;
        }
    } else {
        //VectorXd d_l = this->pfm->electrolyte_diffuse.initial_node->eval(vars);
        for (int i = 0; i < n * mesh.separator_elements.size(); i++) {
            arr_d_eff[i] = pfm->f_diffuse_l(vars(i, 0), vars(i, 1)) / d_ref * eff_mat;
        }
    }

    int i = 0;
    for (auto e: mesh.separator_elements) {
        MatrixXd e_p(n, 1);
        MatrixXd e_c(n, 1);
        MatrixXd e_t(n, 1);
        MatrixXd e_dc(n, 1);
        MatrixXd e_dt(n, 1);

        for (int j = 0; j < n; j++) {
            std::size_t node_id = mesh.elements[e * n + j];
            e_p(j, 0) = u(dof.get_dof(node_id, 0), 0);
            e_c(j, 0) = u(dof.get_dof(node_id, 1), 0);
            e_dc(j, 0) = du(dof.get_dof(node_id, 1), 0);
            if constexpr(use_temp) {
                e_t(j, 0) = u(dof.get_dof(node_id, 4), 0);
                e_dt(j, 0) = du(dof.get_dof(node_id, 4), 0);
            } else {
                e_t(j, 0) = constant::T;
                e_dt(j, 0) = 0;
            }
        }
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

        for (int j = 0; j < n; j++) {
            const MatrixXd &N = shapes.cached_matrix_N[e * n + j];
            const MatrixXd &dN = shapes.cached_matrix_dN[e * n + j];
            //const MatrixXd &NdNT = shapes.cached_matrix_NdNT[e * n + j];
            const MatrixXd &NNT = shapes.cached_matrix_NNT[e * n + j];
            const MatrixXd &dNdNT = shapes.cached_matrix_dNdNT[e * n + j];
            MatrixXd N_T = N.transpose();
            MatrixXd dN_T = dN.transpose();
            double det = shapes.cached_det_J[e * n + j];
            MatrixXd t_mat = N_T * e_c;
            double ele_c_e = t_mat.sum();
            t_mat = N.transpose() * e_t;
            double ele_c_t = t_mat.sum();

            double k_ref = constant::k_ref;
            double k_eff = arr_kappa[i * n + j];
            double dk_dt = 2 * k_eff * constant::R / constant::F * (1 - constant::trans);
            double kd_eff = dk_dt * ele_c_t;
            double d_k_eff = arr_d_kappa[i * n + j];
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
            e_rc += epsilon * eff_1 * NNT * e_dc * w(j) * det + arr_d_eff[i * n + j] * dNdNT * e_c * w(j) * det;
            e_kcc += epsilon * eff_1 * NNT * w(j) * det + arr_d_eff[i * n + j] * dNdNT * w(j) * det;
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
                e_ktt += rho * cap * constant::l_ref * constant::l_ref * NNT / dt * w(j) * det
                         + lambda * dNdNT * w(j) * det;
                e_rt += rho * cap * constant::l_ref * constant::l_ref * NNT * e_dt / dt * w(j) * det
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

        for (int j = 0; j < n; j++) {
            std::size_t id_l = mesh.elements[e * n + j];
            for (int l = 0; l < n; l++) {
                std::size_t id_r = mesh.elements[e * n + l];
                t.emplace_back(dof.get_dof(id_l, 0), dof.get_dof(id_r, 0), e_kpp(j, l));
                t.emplace_back(dof.get_dof(id_l, 0), dof.get_dof(id_r, 1), e_kpc(j, l));
                t.emplace_back(dof.get_dof(id_l, 1), dof.get_dof(id_r, 1), e_kcc(j, l));

                if constexpr(use_temp) {
                    t.emplace_back(dof.get_dof(id_l, 4), dof.get_dof(id_r, 0), e_ktp(j, l));
                    t.emplace_back(dof.get_dof(id_l, 4), dof.get_dof(id_r, 1), e_ktc(j, l));
                    t.emplace_back(dof.get_dof(id_l, 4), dof.get_dof(id_r, 4), e_ktt(j, l));
                }
            }

            res(dof.get_dof(id_l, 0)) += e_rp(j);
            res(dof.get_dof(id_l, 1)) += e_rc(j);
            if constexpr(use_temp) res(dof.get_dof(id_l, 4)) += e_rt(j);
        }
        i++;
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