#include "stiffness_cathode.h"
#include "../../constants/constant.h"
#include "../../functions/functions.h"
#include <eigen3/Eigen/src/Core/Matrix.h>

void stiffness_cathode::generate(const Eigen::Ref<MatrixXd> &u, 
                                 const Eigen::Ref<MatrixXd> &du, 
                                 const Eigen::Ref<MatrixXd> &c_s,
                                 std::vector<Eigen::Triplet<double>> &t, 
                                 Eigen::Ref<VectorXd> res, 
                                 bool is_first_step) {
    const int dim = 1, n = 2;
    int dof_cnt = this->surface_ca_coll - this->surface_an_coll + 1;
    int dof_cnt_eff = dof_cnt - (this->surface_ca_sep - this->surface_an_sep - 1);
    int elem_cnt = dof_cnt - 1;
    int particle_elem_cnt = constant::particle_segment;

    MatrixXd xs = get_integration_point<dim, n>();
    MatrixXd w = get_integration_weight<dim, n>();
    const double dt = constant::dt;
    const double R_p = constant::r_p;
    const double a = 3 * constant::epsilon_s_ca / constant::r_p;
    const double eff_mat = std::pow(constant::epsilon_e_ca, constant::bruggeman);
    const double eff_mat_s = std::pow(constant::epsilon_s_ca, constant::bruggeman);
    const double d_ref = constant::de_ca;
    const double d_eff = constant::de_ca / d_ref * eff_mat;
    const double ds_eff = constant::ds_ca / d_ref;
    const double sigma_ref = constant::sigma_ca * eff_mat_s;
    const double sigma_eff = constant::sigma_ca / sigma_ref * eff_mat_s;
    const double epsilon = constant::epsilon_e_ca;
    const double epsilon_s = constant::epsilon_s_ca;
    const double c_max = constant::c_max_ca;
    const double ce_int = constant::ce_int;
    const double j_ref = constant::j_ref;
    const double F = constant::F;
    const double ce_root = std::sqrt(ce_int);
    const double rho = 2500, cap = 700, lambda = 2.1;

    const int simd_size = this->surface_ca_coll - this->surface_ca_sep + 1;

    double *kqp = new double[simd_size];
    double *kqc = new double[simd_size];
    double *kqs = new double[simd_size];
    //double *kqa = new double[simd_size];
    double *kqq = new double[simd_size];
    double *rq = new double[simd_size];
    //double *rv = new double[simd_size];
    //double *ra = new double[simd_size];
    auto u_ptr = u.data();
    auto du_ptr = du.data();
    std::vector<double> arr_j0(simd_size);
    std::vector<double> arr_d_j0_a(simd_size);
    std::vector<double> arr_d_j0_e(simd_size);
    std::vector<double> arr_uoc(simd_size);
    std::vector<double> arr_d_uoc(simd_size);
    std::vector<double> arr_bv(simd_size);
    std::vector<double> arr_d_bv(simd_size);
    VectorXd c_ss(simd_size);

    for(int i = this->surface_ca_sep - this->surface_an_coll; i <= elem_cnt; ++i) {
        const int idx = i - this->surface_ca_sep + this->surface_an_sep + 1;
        const int src_idx = i - this->surface_ca_sep + this->surface_an_coll;

        c_ss[src_idx] = c_s((idx + 1) * (particle_elem_cnt + 1) - 1, 0);

        arr_j0[src_idx] = j0<2>(u_ptr[dof_cnt + i], c_ss[src_idx]);
        arr_d_j0_a[src_idx] = d_j0_a<2>(u_ptr[dof_cnt + i], c_ss[src_idx]);
        arr_d_j0_e[src_idx] = d_j0_e<2>(u_ptr[dof_cnt + i], c_ss[src_idx]);
        
    }
    if(!settings::use_customize_uoc) {
        for(int i = this->surface_ca_sep - this->surface_an_coll; i <= elem_cnt; ++i) {
            const int idx = i - this->surface_ca_sep + this->surface_an_sep + 1;
            const int src_idx = i - this->surface_ca_sep + this->surface_an_coll;

            arr_uoc[src_idx] = uoc<2>(c_ss[src_idx]);
            arr_d_uoc[src_idx] = d_uoc<2>(c_ss[src_idx]);
            arr_bv[src_idx] = bv(u_ptr[2 * dof_cnt + idx] - u_ptr[i] - arr_uoc[src_idx]);
            arr_d_bv[src_idx] = d_bv(u_ptr[2 * dof_cnt + idx] - u_ptr[i] - arr_uoc[src_idx]);
        }
    } else {
        VectorXd uoc = this->pfm->uoc_cathode.initial_node->eval(c_ss);
        VectorXd d_uoc = this->pfm->uoc_cathode.initial_node->eval_deriv(c_ss, 0);
        for(int i = this->surface_ca_sep - this->surface_an_coll; i <= elem_cnt; ++i) {
            const int idx = i - this->surface_ca_sep + this->surface_an_sep + 1;
            const int src_idx = i - this->surface_ca_sep + this->surface_an_coll;

            arr_uoc[src_idx] = uoc(src_idx);
            arr_d_uoc[src_idx] = d_uoc(src_idx);
            arr_bv[src_idx] = bv(u_ptr[2 * dof_cnt + idx] - u_ptr[i] - arr_uoc[src_idx]);
            arr_d_bv[src_idx] = d_bv(u_ptr[2 * dof_cnt + idx] - u_ptr[i] - arr_uoc[src_idx]);
        }
    }

    for(int i = this->surface_ca_sep - this->surface_an_coll; i <= elem_cnt; ++i) {
        const int idx = i - this->surface_ca_sep + this->surface_an_sep + 1;
        const int src_idx = i - this->surface_ca_sep + this->surface_an_coll;
        double j0_v = arr_j0[src_idx];
        double d_j0_a_v = arr_d_j0_a[src_idx];
        double d_j0_e_v = arr_d_j0_e[src_idx];
        double uoc_v = arr_uoc[src_idx];
        double d_uoc_v = arr_d_uoc[src_idx];
        double bv_v = arr_bv[src_idx];
        double d_bv_v = arr_d_bv[src_idx];

        kqp[src_idx] = -j0_v * d_bv_v * (-1) * c_max * ce_root;
        kqc[src_idx] = -d_j0_e_v * bv_v * c_max * ce_root;
        kqs[src_idx] = -j0_v * d_bv_v * c_max * ce_root;
        kqq[src_idx] = j_ref - (d_j0_a_v * bv_v * c_max * ce_root - j0_v * d_bv_v * d_uoc_v * c_max * ce_root) * dc_ssdj;
        rq[src_idx] = u_ptr[2 * dof_cnt + dof_cnt_eff + idx] * j_ref - j0_v * bv_v * c_max * ce_root;
    }


    for(int i = this->surface_ca_sep - this->surface_an_coll; i <= elem_cnt; ++i) {
        const int idx = i - this->surface_ca_sep + this->surface_an_sep + 1;
        const int src_idx = i - this->surface_ca_sep + this->surface_an_coll;

        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, i, kqp[src_idx]);
        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, i + dof_cnt, kqc[src_idx]);
        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, idx + 2 * dof_cnt, kqs[src_idx]);
        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, idx + 2 * dof_cnt + dof_cnt_eff, kqq[src_idx]);

        res(idx + 2 * dof_cnt + dof_cnt_eff) = rq[src_idx];
    }

    for(int i = this->surface_ca_sep - this->surface_an_coll; i < elem_cnt; i++) {
        int idx = i - this->surface_ca_sep + this->surface_an_sep + 1;
        MatrixXd e_p = u({i, i + 1}, 0); // phi_e
        MatrixXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0); // c_e
        MatrixXd e_s = u({2 * dof_cnt + idx, 2 * dof_cnt + idx + 1}, 0); // phi_s
        MatrixXd e_q = u({2 * dof_cnt + dof_cnt_eff + idx, 2 * dof_cnt + dof_cnt_eff + idx + 1}, 0); // j
        MatrixXd e_t = u({2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll, 2 * dof_cnt + 2 * dof_cnt_eff + i + 1 + this->surface_an_coll}, 0);

        MatrixXd e_dc = du({dof_cnt + i, dof_cnt + i + 1}, 0);
        MatrixXd e_dt = du({2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll, 2 * dof_cnt + 2 * dof_cnt_eff + i + 1 + this->surface_an_coll}, 0);

        MatrixXd e_kss = MatrixXd::Zero(n, n);
        MatrixXd e_ksq = MatrixXd::Zero(n, n);
        MatrixXd e_kcc = MatrixXd::Zero(n, n);
        MatrixXd e_kcq = MatrixXd::Zero(n, n);
        MatrixXd e_kpp = MatrixXd::Zero(n, n);
        MatrixXd e_kpc = MatrixXd::Zero(n, n);
        MatrixXd e_kpq = MatrixXd::Zero(n, n);
        MatrixXd e_ktt = MatrixXd::Zero(n, n);

        MatrixXd e_rs = MatrixXd::Zero(n, 1);
        MatrixXd e_rc = MatrixXd::Zero(n, 1);
        MatrixXd e_rp = MatrixXd::Zero(n, 1);
        MatrixXd e_rt = MatrixXd::Zero(n, 1);

        for(int j = 0; j < n; j++) {
            const MatrixXd &N = cached_matrix_N[(i + surface_an_coll) * n + j];
            const MatrixXd &dN = cached_matrix_dN[(i + surface_an_coll) * n + j];
            const MatrixXd &NNT = cached_matrix_NNT[(i + surface_an_coll) * n + j];
            const MatrixXd &dNdNT = cached_matrix_dNdNT[(i + surface_an_coll) * n + j];
            MatrixXd N_T = N.transpose();
            MatrixXd dN_T = dN.transpose();
            double det = cached_det_J[(i + surface_an_coll) * n + j];
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
            e_kss += sigma_eff * dNdNT * w(j) * det;
            e_ksq += eff_1 * NNT * w(j) * det * j_ref;
            e_rs += sigma_eff * dNdNT * e_s * w(j) * det + eff_1 * NNT * e_q * w(j) * det * j_ref;
            //e_kss = MatrixXd::Identity(2, 2);

            // c part
            double eff_2 = a * constant::l_ref * constant::l_ref * (1 - constant::trans) / d_ref;
            double eff_3 = 1 / dt * constant::l_ref * constant::l_ref / d_ref;
            e_kcc += epsilon * eff_3 * NNT * w(j) * det + d_eff * dNdNT * w(j) * det;
            e_kcq += -eff_2 * NNT * w(j) * det * j_ref / ce_int;
            e_rc += epsilon * eff_3 * NNT * e_dc * w(j) * det + d_eff * dNdNT * e_c * w(j) * det - eff_2 * NNT * e_q * w(j) * det * j_ref / ce_int;
            //e_kcc = MatrixXd::Identity(2, 2);

            // p part
            double eff_4 = a * F * constant::l_ref * constant::l_ref / k_ref;
            e_kpp += k_eff * dNdNT * w(j) * det;
            e_kpc += d_k_eff * dNdNT * e_p * N_T * w(j) * det 
                     - kd_eff / lower * dNdNT * w(j) * det
                     - d_kd_eff / lower * dNdNT * e_c * N_T * w(j) * det
                     + kd_eff / (lower * lower) * dNdNT * e_c * N_T * w(j) * det;
            e_kpq += -eff_4 * NNT * w(j) * det * j_ref;
            e_rp += k_eff * dNdNT * e_p * w(j) * det - kd_eff / lower * dNdNT * e_c * w(j) * det - eff_4 * NNT * e_q * w(j) * det * j_ref;
            //e_kpp = MatrixXd::Identity(2, 2);

            // t part
            e_ktt += rho * cap * NNT / constant::dt * w(j) * det + lambda * dNdNT * w(j) * det;
            e_rt += rho * cap * NNT * e_dt / constant::dt * w(j) * det + lambda * dNdNT * e_t * w(j) * det;
        }

        for(int j = 0; j < n; j++) {
            for(int l = 0; l < n; l++) {
                t.emplace_back(i + j, i + l, e_kpp(j, l));
                t.emplace_back(i + j, i + l + dof_cnt, e_kpc(j, l));
                t.emplace_back(i + j, idx + l + 2 * dof_cnt + dof_cnt_eff, e_kpq(j, l));

                t.emplace_back(i + j + dof_cnt, i + l + dof_cnt, e_kcc(j, l));
                t.emplace_back(i + j + dof_cnt, idx + l + 2 * dof_cnt + dof_cnt_eff, e_kcq(j, l));

                if(true) {
                    t.emplace_back(idx + j + 2 * dof_cnt, idx + l + 2 * dof_cnt, e_kss(j, l));
                    t.emplace_back(idx + j + 2 * dof_cnt, idx + l + 2 * dof_cnt + dof_cnt_eff, e_ksq(j, l));
                }

                t.push_back(Eigen::Triplet<double>(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j, 2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + l, e_ktt(j, l)));
            }
            res(i + j) += e_rp(j);
            res(i + j + dof_cnt) += e_rc(j);
            res(idx + j + 2 * dof_cnt) += e_rs(j);
            res(2 * dof_cnt + 2 * dof_cnt_eff + i + this->surface_an_coll + j) += e_rt(j);
        }

    }

    delete[] kqp;
    delete[] kqc;
    delete[] kqs;
    delete[] kqq;
    delete[] rq;
}