#include "stiffness_anode.h"
#include "../../functions/functions.h"
#include "../../constants/constant.h"
#include <cstdio>

MatrixXd matrix_square(MatrixXd in) {
    return std::move(in.array() * in.array());
}

template<bool use_temp>
void stiffness_anode::generate(const Eigen::Ref<MatrixXd> &u,
                               const Eigen::Ref<MatrixXd> &du,
                               const Eigen::Ref<MatrixXd> &c_s,
                               std::vector<Eigen::Triplet<double> > &t,
                               Eigen::Ref<VectorXd> res,
                               bool is_first_step) {
    const int dim = 1, n = 2;
    //int dof_cnt = this->points.size();
    int dof_cnt = this->surface_ca_coll - this->surface_an_coll + 1;
    int elem_cnt = dof_cnt - 1;
    int dof_cnt_eff = dof_cnt - (this->surface_ca_sep - this->surface_an_sep - 1);
    //int elem_cnt = this->points.size() - 1;
    int particle_elem_cnt = constant::particle_segment;
    MatrixXd xs = get_integration_point<dim, n>();
    MatrixXd w = get_integration_weight<dim, n>();
    const double dt = constant::dt;
    const double R_p = constant::r_p;
    const double a = 3 * constant::epsilon_s_an / constant::r_p;
    const double eff_mat = std::pow(constant::epsilon_e_an, constant::bruggeman);
    const double eff_mat_s = std::pow(constant::epsilon_s_an, constant::bruggeman);
    const double d_ref = constant::de_an;
    //const double d_eff = constant::de_an / d_ref * eff_mat;
    const double ds_eff = constant::ds_an / d_ref;
    const double sigma_ref = constant::sigma_an * eff_mat_s;
    const double sigma_eff = constant::sigma_an / sigma_ref * eff_mat_s;
    const double epsilon = constant::epsilon_e_an;
    const double epsilon_s = constant::epsilon_s_an;
    const double c_max = constant::c_max_an;
    const double ce_int = constant::ce_int;
    const double j_ref = constant::j_ref;
    const double F = constant::F;
    const double ce_root = std::sqrt(ce_int);
    const double rho = constant::density_an, cap = constant::capacity_an, lambda = constant::lambda_an;

    // Here we take out the evaluation of the functions for vectorization purpose.
    // TODO: Extract the code below to one method
    // TODO: Switch all functions to functions accepting VectorXd and use VectorXd
    const int simd_size = this->surface_an_sep - this->surface_an_coll + 1;

    double *kqp = new double[simd_size];
    double *kqc = new double[simd_size];
    double *kqs = new double[simd_size];
    double *kqq = new double[simd_size];
    double *rq = new double[simd_size];
    auto u_ptr = u.data();
    auto du_ptr = du.data();
    std::vector<double> arr_j0(simd_size);
    std::vector<double> arr_d_j0_a(simd_size);
    std::vector<double> arr_d_j0_e(simd_size);
    std::vector<double> arr_uoc(simd_size);
    std::vector<double> arr_d_uoc(simd_size);
    std::vector<double> arr_bv(simd_size);
    std::vector<double> arr_d_bv(simd_size);
    std::vector<double> arr_kappa(2 * (simd_size - 1));
    std::vector<double> arr_d_kappa(2 * (simd_size - 1));
    std::vector<double> arr_d_eff(2 * (simd_size - 1));
    double *arr_eta = new double[simd_size];
    double *arr_du = new double[simd_size];
    VectorXd c_ss(simd_size);


    for (int i = 0; i <= this->surface_an_sep - this->surface_an_coll; ++i) {
        const int idx = i;

        c_ss[i] = c_s((idx + 1) * (particle_elem_cnt + 1) - 1, 0);

        arr_j0[i] = j0<1>(u_ptr[dof_cnt + idx], c_ss[i]);
        arr_d_j0_a[i] = d_j0_a<1>(u_ptr[dof_cnt + idx], c_ss[i]);
        arr_d_j0_e[i] = d_j0_e<1>(u_ptr[dof_cnt + idx], c_ss[i]);
    }
    if (!settings::use_customize_uoc) {
        for (int i = 0; i <= this->surface_an_sep - this->surface_an_coll; ++i) {
            const int idx = i;

            arr_uoc[i] = uoc<1>(c_ss[i]);
            arr_eta[i] = u_ptr[2 * dof_cnt + idx] - u_ptr[i] - arr_uoc[i];
            arr_d_uoc[i] = d_uoc<1>(c_ss[i]);
            arr_bv[i] = bv(arr_eta[i]);
            arr_d_bv[i] = d_bv(arr_eta[i]);
            arr_du[i] = arr_d_uoc[i] / c_max;
        }
    } else {
        VectorXd uoc = this->pfm->uoc_anode.initial_node->eval(c_ss);
        VectorXd d_uoc = this->pfm->uoc_anode.initial_node->eval_deriv(c_ss, 0);
        VectorXd v_du = this->pfm->anode_entropy.initial_node->eval(c_ss);
        for (int i = 0; i <= this->surface_an_sep - this->surface_an_coll; ++i) {
            const int idx = i;

            arr_uoc[i] = uoc(i);
            arr_eta[i] = u_ptr[2 * dof_cnt + idx] - u_ptr[i] - arr_uoc[i];
            arr_d_uoc[i] = d_uoc(i);
            arr_bv[i] = bv(arr_eta[i]);
            arr_d_bv[i] = d_bv(arr_eta[i]);
            arr_du[i] = v_du(i);
        }
    }

    MatrixXd vars(2 * (simd_size - 1), 2);
    if constexpr (use_temp) {
        for (int i = 0; i < this->surface_an_sep - this->surface_an_coll; ++i) {
            VectorXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0);
            VectorXd e_t = u({
                                 2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll,
                                 2 * dof_cnt + 2 * dof_cnt_eff + i + 1 + this->surface_an_coll
                             }, 0);
            for (int j = 0; j < n; j++) {
                const MatrixXd &N = cached_matrix_N[(i + surface_an_coll) * n + j];

                MatrixXd t_mat = N.transpose() * e_c;
                double ele_c_e = t_mat.sum();
                t_mat = N.transpose() * e_t;
                double ele_c_t = t_mat.sum();
                vars(i * n + j, 0) = ele_c_e * ce_int;
                vars(i * n + j, 1) = ele_c_t;
            }
        }
    } else {
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
    }
    if (!settings::use_customize_kappa) {
        for (int i = 0; i < this->surface_an_sep - this->surface_an_coll; ++i) {
            for (int j = 0; j < n; j++) {
                double k_eff = kappa(vars(i * n + j, 0)) / constant::k_ref * eff_mat;
                double d_k_eff = d_kappa(vars(i * n + j, 0)) * ce_int / constant::k_ref * eff_mat;
                arr_kappa[i * n + j] = k_eff;
                arr_d_kappa[i * n + j] = d_k_eff;
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
        for (int i = 0; i < this->surface_an_sep - this->surface_an_coll; ++i) {
            for (int j = 0; j < n; j++) {
                arr_d_eff[i * n + j] = constant::de_an / d_ref * eff_mat;
            }
        }
    } else {
        for (int i = 0; i < this->surface_an_sep - this->surface_an_coll; ++i) {
            VectorXd d_l = this->pfm->electrolyte_diffuse.initial_node->eval(vars);
            for (int j = 0; j < n; j++) {
                arr_d_eff[i * n + j] = d_l(i * n + j) / d_ref * eff_mat;
            }
        }
    }

    for (int i = 0; i <= this->surface_an_sep - this->surface_an_coll; ++i) {
        const int idx = i;

        double j0_v = arr_j0[idx];
        double d_j0_a_v = arr_d_j0_a[idx];
        double d_j0_e_v = arr_d_j0_e[idx];
        double uoc_v = arr_uoc[idx];
        double d_uoc_v = arr_d_uoc[idx];
        double bv_v = arr_bv[idx];
        double d_bv_v = arr_d_bv[idx];

        kqp[i] = -j0_v * d_bv_v * (-1) * c_max * ce_root;
        kqc[i] = -d_j0_e_v * bv_v * c_max * ce_root;
        kqs[i] = -j0_v * d_bv_v * c_max * ce_root;
        kqq[i] = j_ref - (d_j0_a_v * bv_v * c_max * ce_root - j0_v * d_bv_v * d_uoc_v * c_max * ce_root) * dc_ssdj;
        rq[i] = u_ptr[2 * dof_cnt + dof_cnt_eff + idx] * j_ref - j0_v * bv_v * c_max * ce_root;
    }

    for (int i = 0; i <= this->surface_an_sep - this->surface_an_coll; ++i) {
        const int idx = i;
        const int src_idx = i;

        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, i, kqp[src_idx]);
        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, i + dof_cnt, kqc[src_idx]);
        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, idx + 2 * dof_cnt, kqs[src_idx]);
        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, idx + 2 * dof_cnt + dof_cnt_eff, kqq[src_idx]);

        res(idx + 2 * dof_cnt + dof_cnt_eff) = rq[src_idx];
    }

    // Total assembly
    for (int i = 0; i < this->surface_an_sep - this->surface_an_coll; i++) {
        MatrixXd e_p = u({i, i + 1}, 0); // phi_e
        MatrixXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0); // c_e
        MatrixXd e_s = u({2 * dof_cnt + i, 2 * dof_cnt + i + 1}, 0); // phi_s
        MatrixXd e_q = u({2 * dof_cnt + dof_cnt_eff + i, 2 * dof_cnt + dof_cnt_eff + i + 1}, 0); // j
        MatrixXd e_t = MatrixXd::Ones(n, 1) * constant::T;
        if constexpr (use_temp) {
            e_t = u({
                                 2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll,
                                 2 * dof_cnt + 2 * dof_cnt_eff + i + 1 + this->surface_an_coll
                             }, 0);
        }

        MatrixXd e_dc = du({dof_cnt + i, dof_cnt + i + 1}, 0);
        MatrixXd e_dt = MatrixXd::Zero(n, 1);
        if constexpr (use_temp) {
            e_dt = du({
                                   2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll,
                                   2 * dof_cnt + 2 * dof_cnt_eff + i + 1 + this->surface_an_coll
                               }, 0);
        }

        MatrixXd e_kss = MatrixXd::Zero(n, n);
        MatrixXd e_ksq = MatrixXd::Zero(n, n);
        MatrixXd e_kcc = MatrixXd::Zero(n, n);
        MatrixXd e_kcq = MatrixXd::Zero(n, n);
        MatrixXd e_kpp = MatrixXd::Zero(n, n);
        MatrixXd e_kpc = MatrixXd::Zero(n, n);
        MatrixXd e_kpq = MatrixXd::Zero(n, n);
        MatrixXd e_ktt = MatrixXd::Zero(n, n);
        MatrixXd e_ktp = MatrixXd::Zero(n, n);
        MatrixXd e_ktc = MatrixXd::Zero(n, n);
        MatrixXd e_ktq = MatrixXd::Zero(n, n);

        MatrixXd e_rs = MatrixXd::Zero(n, 1);
        MatrixXd e_rc = MatrixXd::Zero(n, 1);
        MatrixXd e_rp = MatrixXd::Zero(n, 1);
        MatrixXd e_rt = MatrixXd::Zero(n, 1);

        for (int j = 0; j < n; j++) {
            const MatrixXd &N = cached_matrix_N[(i + surface_an_coll) * n + j];
            const MatrixXd &dN = cached_matrix_dN[(i + surface_an_coll) * n + j];
            const MatrixXd &NdNT = cached_matrix_NdNT[(i + surface_an_coll) * n + j];
            const MatrixXd &NNT = cached_matrix_NNT[(i + surface_an_coll) * n + j];
            const MatrixXd &dNdNT = cached_matrix_dNdNT[(i + surface_an_coll) * n + j];
            MatrixXd N_T = N.transpose();
            MatrixXd dN_T = dN.transpose();
            double det = cached_det_J[(i + surface_an_coll) * n + j];
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
            double d_eff = arr_d_eff[i * n + j];

            // s part
            double eff_1 = a * F * constant::l_ref * constant::l_ref / sigma_ref;
            e_kss += sigma_eff * dNdNT * w(j) * det;
            e_ksq += eff_1 * NNT * w(j) * det * j_ref;
            e_rs += sigma_eff * dNdNT * e_s * w(j) * det + eff_1 * NNT * e_q * w(j) * det * j_ref;
            //e_kss = MatrixXd::Identity(2, 2);

            // c part
            double eff_2 = a * constant::l_ref * constant::l_ref * (1 - constant::trans) / d_ref;
            double eff_3 = 1 / dt * constant::l_ref * constant::l_ref / d_ref;
            e_kcc += epsilon * eff_3 * NNT * w(j) * det + d_eff * dNdNT * w(j) * det;
            e_kcq += -eff_2 * NNT * w(j) * det * j_ref / ce_int;
            e_rc += epsilon * eff_3 * NNT * e_dc * w(j) * det + d_eff * dNdNT * e_c * w(j) * det - eff_2 * NNT * e_q *
                    w(j) * det * j_ref / ce_int;
            //e_kcc = MatrixXd::Identity(2, 2);

            // p part
            double eff_4 = a * F * constant::l_ref * constant::l_ref / k_ref;
            e_kpp += k_eff * dNdNT * w(j) * det;
            e_kpc += d_k_eff * dNdNT * e_p * N_T * w(j) * det
                    - kd_eff / ele_c_e * dNdNT * w(j) * det
                    - d_kd_eff / ele_c_e * dNdNT * e_c * N_T * w(j) * det
                    + kd_eff / (ele_c_e * ele_c_e) * dNdNT * e_c * N_T * w(j) * det;
            e_kpq += -eff_4 * NNT * w(j) * det * j_ref;
            e_rp += k_eff * dNdNT * e_p * w(j) * det - kd_eff / ele_c_e * dNdNT * e_c * w(j) * det - eff_4 * NNT * e_q *
                    w(j) * det * j_ref;
            //e_kpp = MatrixXd::Identity(2, 2);

            // t part
            if constexpr (use_temp) {
                MatrixXd e_eta = Eigen::Map<VectorXd>(arr_eta + i, 2);
                MatrixXd e_du = Eigen::Map<VectorXd>(arr_du + i, 2);
                double dNe_p = (dN_T * e_p).sum();
                double dNe_s = (dN_T * e_s).sum();
                double Ne_t = (N_T * e_t).sum();
                double dNe_c = (dN_T * e_c).sum();
                double Ne_q = (N_T * e_q).sum();
                double Ne_eta = (N_T * e_eta).sum();
                double Ne_du = (N_T * e_du).sum();
                double e_dp2 = dNe_p * dNe_p;
                double e_ds2 = dNe_s * dNe_s;
                double e_tdpdc = Ne_t * dNe_p * dNe_c;
                double e_dpdc = dNe_p * dNe_c;
                double e_tdc = Ne_t * dNe_c;
                double e_tdp = Ne_t * dNe_p;
                double e_qt = Ne_q * Ne_t;
                double e_qeta = Ne_q * Ne_eta;
                double e_qtdu = Ne_q * Ne_t * Ne_du;
                double e_qdu = Ne_q * Ne_du;
                // Heat transfer
                e_ktt += rho * cap * constant::l_ref * constant::l_ref * NNT / constant::dt * w(j) * det
                        + lambda * dNdNT * w(j) * det;
                e_rt += rho * cap * constant::l_ref * constant::l_ref * NNT * e_dt / constant::dt * w(j) * det
                        + lambda * dNdNT * e_t * w(j) * det;
                // Q_ohm
                e_ktt += -(dk_dt * NNT * e_dpdc / ele_c_e);
                e_ktp += -(k_eff * NdNT * 2 * dNe_p + dk_dt * NdNT * e_tdc / ele_c_e);
                e_ktc += -(dk_dt * NdNT * e_tdp / ele_c_e);
                e_rt += -(k_eff * N * e_dp2 + dk_dt * N * e_tdpdc / ele_c_e) * k_ref * w(j) * det;
                e_rt += -(sigma_eff * N * e_ds2) * sigma_ref * w(j) * det;
                // Q_rxn
                e_ktq += -(constant::F * a * NNT * Ne_eta) * constant::l_ref * constant::l_ref * j_ref * w(j) * det;
                e_rt += -(constant::F * a * N * e_qeta) * constant::l_ref * constant::l_ref * j_ref * w(j) * det;
                // Q_rev
                e_ktt += -(constant::F * a * NNT * e_qdu) * constant::l_ref * constant::l_ref * j_ref * w(j) * det;
                e_rt += -(constant::F * a * N * e_qtdu) * constant::l_ref * constant::l_ref * j_ref * w(j) * det;
            }
        }

        for (int j = 0; j < n; j++) {
            for (int l = 0; l < n; l++) {
                if (true) {
                    t.emplace_back(i + j, i + l, e_kpp(j, l));
                    t.emplace_back(i + j, i + l + dof_cnt, e_kpc(j, l));
                    t.emplace_back(i + j, i + l + 2 * dof_cnt + dof_cnt_eff, e_kpq(j, l));
                }

                t.emplace_back(i + j + dof_cnt, i + l + dof_cnt, e_kcc(j, l));
                t.emplace_back(i + j + dof_cnt, i + l + 2 * dof_cnt + dof_cnt_eff, e_kcq(j, l));

                if (i != 0 || j != 0) {
                    t.emplace_back(i + j + 2 * dof_cnt, i + l + 2 * dof_cnt, e_kss(j, l));
                    t.emplace_back(i + j + 2 * dof_cnt, i + l + 2 * dof_cnt + dof_cnt_eff, e_ksq(j, l));
                }

                if constexpr(use_temp) {
                    t.emplace_back(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j,
                                   2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + l, e_ktt(j, l));
                    t.emplace_back(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j, i + l, e_ktp(j, l));
                    t.emplace_back(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j, i + l + dof_cnt,
                                   e_ktc(j, l));
                    t.emplace_back(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j,
                                   i + l + 2 * dof_cnt + dof_cnt_eff, e_ktq(j, l));
                }
            }
            res(i + j) += e_rp(j);
            res(i + j + dof_cnt) += e_rc(j);
            res(i + j + 2 * dof_cnt) += e_rs(j);
            if constexpr (use_temp) res(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j) += e_rt(j);
        }
    }

    delete[] kqp;
    delete[] kqc;
    delete[] kqs;
    delete[] kqq;
    delete[] rq;
    delete[] arr_eta;
    delete[] arr_du;
}

template void stiffness_anode::generate<true>(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                                              const Eigen::Ref<MatrixXd> &,
                                              std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                                              bool);

template void stiffness_anode::generate<false>(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                                               const Eigen::Ref<MatrixXd> &,
                                               std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                                               bool);
