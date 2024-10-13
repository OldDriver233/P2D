#include "stiffness_anode.h"
#include "../../functions/functions.h"
#include "../../constants/constant.h"

void stiffness_anode::generate(const Eigen::Ref<MatrixXd> &u, const Eigen::Ref<MatrixXd> &du, std::vector<Eigen::Triplet<double>> &t, Eigen::Ref<VectorXd> res, bool is_first_step) {
    const int dim = 1, n = 2;
    int dof_cnt = this->points.size();
    int dof_cnt_eff = dof_cnt - (this->surface_ca_sep - this->surface_an_sep - 1);
    int elem_cnt = this->points.size() - 1;
    MatrixXd xs = get_integration_point<dim, n>();
    MatrixXd w = get_integration_weight<dim, n>();
    const double dt = constant::dt;
    const double R_p = constant::r_p;
    const double a = 3 * constant::epsilon_s_an / constant::r_p;
    const double eff_mat = std::pow(constant::epsilon_e_an, constant::bruggeman);
    const double eff_mat_s = std::pow(constant::epsilon_s_an, constant::bruggeman);
    const double d_ref = constant::de_an;
    const double d_eff = constant::de_an / d_ref * eff_mat;
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

    const int simd_size = this->surface_an_sep + 1;

    double *kqp = new double[simd_size];
    double *kqc = new double[simd_size];
    double *kqs = new double[simd_size];
    double *kqa = new double[simd_size];
    double *rq = new double[simd_size];
    double *rv = new double[simd_size];
    double *ra = new double[simd_size];
    auto u_ptr = u.data();
    auto du_ptr = du.data();

    #pragma omp simd
    for(int i = 0; i <= this->surface_an_sep; ++i) {
        const int idx = i;

        const double e_dv = du_ptr[2 * dof_cnt + 2 * dof_cnt_eff + idx];
        const double e_q = u_ptr[2 * dof_cnt + dof_cnt_eff + idx];
        const double e_v = u_ptr[2 * dof_cnt + 2 * dof_cnt_eff + idx];
        const double e_a = u_ptr[2 * dof_cnt + 3 * dof_cnt_eff + idx];

        double j0_v = j0<1>(u_ptr[dof_cnt + idx], u_ptr[2 * dof_cnt + 3 * dof_cnt_eff + idx]);
        double d_j0_a_v = d_j0_a<1>(u_ptr[dof_cnt + idx], u_ptr[2 * dof_cnt + 3 * dof_cnt_eff + idx]);
        double d_j0_e_v = d_j0_e<1>(u_ptr[dof_cnt + idx], u_ptr[2 * dof_cnt + 3 * dof_cnt_eff + idx]);
        double uoc_v = uoc<1>(u_ptr[2 * dof_cnt + 3 * dof_cnt_eff + idx]);
        double d_uoc_v = d_uoc<1>(u_ptr[2 * dof_cnt + 3 * dof_cnt_eff + idx]);
        double bv_v = bv(u_ptr[2 * dof_cnt + idx] - u_ptr[i] - uoc_v);
        double d_bv_v = d_bv(u_ptr[2 * dof_cnt + idx] - u_ptr[i] - uoc_v);

        kqp[i] = -j0_v * d_bv_v * (-1) * c_max * ce_root / j_ref;
        kqc[i] = -d_j0_e_v * bv_v * c_max * ce_root / j_ref;
        kqs[i] = -j0_v * d_bv_v * c_max * ce_root / j_ref;
        kqa[i] = -d_j0_a_v * bv_v * c_max * ce_root / j_ref - j0_v * d_bv_v * (-1) * d_uoc_v * c_max * ce_root / j_ref;
        rq[i] = u_ptr[2 * dof_cnt + dof_cnt_eff + idx] - j0_v * bv_v * c_max * ce_root / j_ref;
        rv[i] = e_dv / dt + 3 * e_q / R_p / c_max * j_ref;
        ra[i] = e_a - e_v + e_q * 0.2 * R_p / ds_eff / c_max * j_ref / d_ref;
    }

    const double kvq = 3 / R_p / c_max * j_ref;
    const double kvv = 1 / dt;
    const double kaa = 1;
    const double kav = -1;
    const double kaq = 0.2 * R_p / ds_eff / c_max * j_ref / d_ref;
    const double kqq = 1;

    for(int i = 0; i <= this->surface_an_sep; ++i) {
        const int idx = i;
        const int src_idx = i;

        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, i, kqp[src_idx]);
        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, i + dof_cnt, kqc[src_idx]);
        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, idx + 2 * dof_cnt, kqs[src_idx]);
        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, idx + 2 * dof_cnt + dof_cnt_eff, kqq);
        t.emplace_back(2 * dof_cnt + dof_cnt_eff + idx, idx + 2 * dof_cnt + 3 * dof_cnt_eff, kqa[src_idx]);

        t.emplace_back(2 * dof_cnt + 2 * dof_cnt_eff + idx, idx + 2 * dof_cnt + dof_cnt_eff, kvq);
        t.emplace_back(2 * dof_cnt + 2 * dof_cnt_eff + idx, idx + 2 * dof_cnt + 2 * dof_cnt_eff, kvv);

        t.emplace_back(2 * dof_cnt + 3 * dof_cnt_eff + idx, idx + 2 * dof_cnt + dof_cnt_eff, kaq);
        t.emplace_back(2 * dof_cnt + 3 * dof_cnt_eff + idx, idx + 2 * dof_cnt + 2 * dof_cnt_eff, kav);
        t.emplace_back(2 * dof_cnt + 3 * dof_cnt_eff + idx, idx + 2 * dof_cnt + 3 * dof_cnt_eff, kaa);

        res(idx + 2 * dof_cnt + dof_cnt_eff) = rq[src_idx];
        res(idx + 2 * dof_cnt + 2 * dof_cnt_eff) = rv[src_idx];
        res(idx + 2 * dof_cnt + 3 * dof_cnt_eff) = ra[src_idx];
    }

    for(int i = 0; i < this->surface_an_sep; i++) {
        MatrixXd e_p = u({i, i + 1}, 0); // phi_e
        MatrixXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0); // c_e
        MatrixXd e_s = u({2 * dof_cnt + i, 2 * dof_cnt + i + 1}, 0); // phi_s
        MatrixXd e_q = u({2 * dof_cnt + dof_cnt_eff + i, 2 * dof_cnt + dof_cnt_eff + i + 1}, 0); // j
        MatrixXd e_v = u({2 * dof_cnt + 2 * dof_cnt_eff + i, 2 * dof_cnt + 2 * dof_cnt_eff + i + 1}, 0); // c^avg
        MatrixXd e_a = u({2 * dof_cnt + 3 * dof_cnt_eff + i, 2 * dof_cnt + 3 * dof_cnt_eff + i + 1}, 0); // c^*

        MatrixXd e_dc = du({dof_cnt + i, dof_cnt + i + 1}, 0);
        MatrixXd e_dv = du({2 * dof_cnt + 2 * dof_cnt_eff + i, 2 * dof_cnt + 2 * dof_cnt_eff + i + 1}, 0);

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

        MatrixXd e_rs = MatrixXd::Zero(n, 1);
        MatrixXd e_rc = MatrixXd::Zero(n, 1);
        MatrixXd e_rp = MatrixXd::Zero(n, 1);

        for(int j = 0; j < n; j++) {
            const MatrixXd &N = cached_matrix_N[i * n + j];
            const MatrixXd &dN = cached_matrix_dN[i * n + j];
            const MatrixXd &NNT = cached_matrix_NNT[i * n + j];
            const MatrixXd &dNdNT = cached_matrix_dNdNT[i * n + j];
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
        }

        for(int j = 0; j < n; j++) {
            for(int l = 0; l < n; l++) {
                if((!is_first_step) || i != 0 || j != 0) {
                    t.push_back(Eigen::Triplet<double>(i + j, i + l, e_kpp(j, l)));
                    t.push_back(Eigen::Triplet<double>(i + j, i + l + dof_cnt, e_kpc(j, l)));
                    t.push_back(Eigen::Triplet<double>(i + j, i + l + 2 * dof_cnt + dof_cnt_eff, e_kpq(j, l)));
                }

                t.push_back(Eigen::Triplet<double>(i + j + dof_cnt, i + l + dof_cnt, e_kcc(j, l)));
                t.push_back(Eigen::Triplet<double>(i + j + dof_cnt, i + l + 2 * dof_cnt + dof_cnt_eff, e_kcq(j, l)));

                if(i != 0 || j != 0) {
                    t.push_back(Eigen::Triplet<double>(i + j + 2 * dof_cnt, i + l + 2 * dof_cnt, e_kss(j, l)));
                    t.push_back(Eigen::Triplet<double>(i + j + 2 * dof_cnt, i + l + 2 * dof_cnt + dof_cnt_eff, e_ksq(j, l)));
                }
            }
            res(i + j) += e_rp(j);
            res(i + j + dof_cnt) += e_rc(j);
            res(i + j + 2 * dof_cnt) += e_rs(j);
        }

    }

    delete[] kqp;
    delete[] kqc;
    delete[] kqs;
    delete[] kqa;
    delete[] rq;
    delete[] rv;
    delete[] ra;
}