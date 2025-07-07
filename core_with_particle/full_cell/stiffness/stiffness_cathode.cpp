#include "stiffness_cathode.h"
#include "../../constants/constant.h"
#include "../../functions/functions.h"
#include <eigen3/Eigen/src/Core/Matrix.h>

template<bool use_temp>
void stiffness_cathode::generate(const Eigen::Ref<MatrixXd> &u, 
                                 const Eigen::Ref<MatrixXd> &du, 
                                 const Eigen::Ref<MatrixXd> &c_s,
                                 const std::vector<double> &stress,
                                 const std::vector<double> &avg_conc,
                                 std::vector<Eigen::Triplet<double>> &t, 
                                 Eigen::Ref<VectorXd> res, 
                                 bool is_first_step,
                                 std::vector<double>& temp) {
    const int dim = get_dim(mesh.p_type), n = get_nodes(mesh.p_type);
    std::size_t node_cnt = mesh.node_count;
    int particle_elem_cnt = constant::particle_segment;

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

    const double dt = this->st->dt_now;
    const double R_p = constant::cathode.r_p;
    const double a = 3 * constant::cathode.epsilon_s / R_p;
    const double eff_mat = std::pow(constant::cathode.epsilon_e, constant::bruggeman);
    const double eff_mat_s = std::pow(constant::cathode.epsilon_s, constant::bruggeman);
    const double d_ref = constant::cathode.D_e;
    //const double d_eff = constant::de_ca / d_ref * eff_mat;
    //const double ds_eff = constant::ds_ca / d_ref;
    const double sigma_ref = constant::cathode.sigma * eff_mat_s;
    const double sigma_eff = constant::cathode.sigma / sigma_ref * eff_mat_s;
    const double epsilon = constant::cathode.epsilon_e;
    const double epsilon_s = constant::cathode.epsilon_s;
    const double c_max = constant::cathode.c_max;
    const double ce_int = constant::ce_int;
    const double j_ref = constant::j_ref;
    const double F = constant::F;
    const double ce_root = std::sqrt(ce_int);
    const double rho = constant::cathode.density, cap = constant::cathode.capacity, lambda = constant::cathode.lambda;

    const int simd_size = mesh.cathode_nodes.size();

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
    std::vector<double> arr_kappa(n * mesh.cathode_elements.size());
    std::vector<double> arr_d_kappa(n * mesh.cathode_elements.size());
    std::vector<double> arr_d_eff(n * mesh.cathode_elements.size());
    std::vector<double> ktt(simd_size);
    std::vector<double> ktq(simd_size);
    std::vector<double> rt(simd_size);
    std::vector<double> stress_surf(simd_size);
    double *arr_eta = new double[simd_size];
    double *arr_du = new double[simd_size];
    VectorXd c_ss(simd_size);

    int i = 0;
    for (auto x: mesh.cathode_nodes) {
        double t = constant::T;
        double c_e = 1;
        if (settings::calc_temperature) {
            t = u_ptr[dof.get_dof(x, 4)];
        }
        if (!settings::is_solid_battery) {
            c_e = u_ptr[dof.get_dof(x, 1)];
        }
        c_ss[i] = c_s((dof.particle_mapper[x] + 1) * (particle_elem_cnt + 1) - 1, 0);

        arr_j0[i] = j0<2>(c_e, c_ss[i], t);
        arr_d_j0_a[i] = d_j0_a<2>(c_e, c_ss[i], t);
        arr_d_j0_e[i] = d_j0_e<2>(c_e, c_ss[i], t);
        i++;
    }
    if (settings::stress_analysis) {
        i = 0;
        const double E_p = constant::cathode.E, nu_p = constant::cathode.nu, omega = constant::cathode.omega;
        for (auto x: mesh.cathode_nodes) {
            int i_avg = dof.particle_mapper[x];
            double stress_outer = (stress[4 * x] + stress[4 * x + 1] + stress[4 * x + 2]) / (3 * epsilon_s);
            stress_surf[i] = 2 * omega * E_p / (9 * (1 - nu_p)) * (avg_conc[i_avg] - c_ss[i] * c_max) + stress_outer;
            i++;
        }
    } else {
        i = 0;
        for (auto x: mesh.cathode_nodes) {
            stress_surf[i] = 0;
            i++;
        }
    }
    if (!settings::use_customize_uoc) {
        int i = 0;
        for (auto x: mesh.cathode_nodes) {
            double t = constant::T;
            double omega = constant::cathode.omega;
            if (settings::calc_temperature) {
                t = u_ptr[dof.get_dof(x, 4)];
            }

            arr_uoc[i] = uoc<2>(c_ss[i]);
            arr_eta[i] = u_ptr[dof.get_dof(x, 2)] - u_ptr[dof.get_dof(x, 0)] - arr_uoc[i] - omega * stress_surf[i] / constant::F;
            arr_d_uoc[i] = d_uoc<2>(c_ss[i]);
            arr_bv[i] = bv(arr_eta[i], t);
            arr_d_bv[i] = d_bv(arr_eta[i], t);
            arr_du[i] = arr_d_uoc[i] / c_max;
            i++;
        }
    } else {
        int i = 0;
        for (auto x: mesh.cathode_nodes) {
            double t = constant::T;
            double omega = constant::cathode.omega;

            if (settings::calc_temperature) {
                t = u_ptr[dof.get_dof(x, 4)];
            }
            arr_uoc[i] = this->pfm->f_uoc_cathode(c_ss[i]);
            arr_eta[i] = u_ptr[dof.get_dof(x, 2)] - u_ptr[dof.get_dof(x, 0)] - arr_uoc[i];
            arr_d_uoc[i] = this->pfm->f_d_uoc_cathode(c_ss[i]);
            arr_bv[i] = bv(arr_eta[i], t);
            arr_d_bv[i] = d_bv(arr_eta[i], t);
            arr_du[i] = this->pfm->f_cathode_entropy(c_ss[i]);
            i++;
        }
    }

    MatrixXd vars(n * mesh.cathode_elements.size(), 2);
    if constexpr (use_temp) {
        int i = 0;
        for (auto e: mesh.cathode_elements) {
            VectorXd e_c = VectorXd::Zero(n);
            VectorXd e_t(n);
            for (int j = 0; j < n; j++) {
                int node_id = mesh.elements[e * n + j];
                if (!settings::is_solid_battery) e_c(j) = u(dof.get_dof(node_id, 1), 0);
                e_t(j) = u(dof.get_dof(node_id, 4), 0);
            }
            for (int j = 0; j < n; j++) {
                const MatrixXd &N = shapes.cached_matrix_N[e * n + j];
                MatrixXd t_mat = N.transpose() * e_c;
                double ele_c_e = t_mat.value();
                t_mat = N.transpose() * e_t;
                double ele_c_t = t_mat.value();
                vars(i * n + j, 0) = ele_c_e * ce_int;
                vars(i * n + j, 1) = ele_c_t;
            }
            i++;
        }
    } else {
        int i = 0;
        for (auto e: mesh.cathode_elements) {
            VectorXd e_c = VectorXd::Zero(n);
            for (int j = 0; j < n; j++) {
                int node_id = mesh.elements[e * n + j];
                if (!settings::is_solid_battery) e_c(j) = u(dof.get_dof(node_id, 1), 0);
            }
            for (int j = 0; j < n; j++) {
                const MatrixXd &N = shapes.cached_matrix_N[e * n + j];
                MatrixXd t_mat = N.transpose() * e_c;
                double ele_c_e = t_mat.value();
                vars(i * n + j, 0) = ele_c_e * ce_int;
                vars(i * n + j, 1) = constant::T;
            }
            i++;
        }
    }
    if (!settings::use_customize_kappa) {
        for (int i = 0; i < mesh.cathode_elements.size(); ++i) {
            for (int j = 0; j < n; j++) {
                double k_eff = kappa(vars(i * n + j, 0)) / constant::k_ref * eff_mat;
                double d_k_eff = d_kappa(vars(i * n + j, 0)) * ce_int / constant::k_ref * eff_mat;
                arr_kappa[i * n + j] = k_eff;
                arr_d_kappa[i * n + j] = d_k_eff;
            }
        }
    } else {
        for (int i = 0; i < mesh.cathode_elements.size(); ++i) {
            for (int j = 0; j < n; j++) {
                double k_eff = pfm->f_kappa(vars(i * n + j, 0), vars(i * n + j, 1)) / constant::k_ref * eff_mat;
                double d_k_eff = pfm->f_d_kappa(vars(i * n + j, 0), vars(i * n + j, 1)) * ce_int / constant::k_ref * eff_mat;
                arr_kappa[i * n + j] = k_eff;
                arr_d_kappa[i * n + j] = d_k_eff;
            }
        }
    }
    if (!settings::is_solid_battery) {
        if (!settings::use_customize_diffuse) {
            for (int i = 0; i < mesh.cathode_elements.size(); ++i) {
                for (int j = 0; j < n; j++) {
                    arr_d_eff[i * n + j] = constant::cathode.D_e / d_ref * eff_mat;
                }
            }
        } else {
            for (int i = 0; i < mesh.cathode_elements.size(); ++i) {
                for (int j = 0; j < n; j++) {
                    arr_d_eff[i * n + j] = pfm->f_diffuse_l(vars(i * n + j, 0), vars(i * n + j, 1)) / d_ref * eff_mat;
                }
            }
        }
    }

    i = 0;
    for (auto x: mesh.cathode_nodes) {
        double j0_v = arr_j0[i];
        double d_j0_a_v = arr_d_j0_a[i];
        double d_j0_e_v = arr_d_j0_e[i];
        double uoc_v = arr_uoc[i];
        double d_uoc_v = arr_d_uoc[i];
        double bv_v = arr_bv[i];
        double d_bv_v = arr_d_bv[i];

        if (!settings::is_solid_battery) {
            kqp[i] = -j0_v * d_bv_v * (-1) * c_max * ce_root;
            kqc[i] = -d_j0_e_v * bv_v * c_max * ce_root;
            kqs[i] = -j0_v * d_bv_v * c_max * ce_root;
            kqq[i] = j_ref - (d_j0_a_v * bv_v * c_max * ce_root - j0_v * d_bv_v * d_uoc_v * c_max * ce_root) * dc_ssdj;
            rq[i] = u_ptr[dof.get_dof(x, 3)] * j_ref - j0_v * bv_v * c_max * ce_root;
        } else {
            kqp[i] = -j0_v * d_bv_v * (-1) * c_max;
            kqs[i] = -j0_v * d_bv_v * c_max;
            kqq[i] = j_ref - (d_j0_a_v * bv_v * c_max - j0_v * d_bv_v * d_uoc_v * c_max) * dc_ssdj;
            rq[i] = u_ptr[dof.get_dof(x, 3)] * j_ref - j0_v * bv_v * c_max;
            //std::cout<<u_ptr[dof.get_dof(x, 2)]<<std::endl;
            //std::cout<<u_ptr[dof.get_dof(x, 0)]<<std::endl;
            //std::cout<<rq[i]<<std::endl;
            //std::cout<<bv_v<<std::endl;
            //std::cout<<std::endl;
        }
        i++;
    }


    i = 0;
    for (auto x: mesh.cathode_nodes) {
        t.emplace_back(dof.get_dof(x, 3), dof.get_dof(x, 0), kqp[i]);
        if (!settings::is_solid_battery) t.emplace_back(dof.get_dof(x, 3), dof.get_dof(x, 1), kqc[i]);
        t.emplace_back(dof.get_dof(x, 3), dof.get_dof(x, 2), kqs[i]);
        t.emplace_back(dof.get_dof(x, 3), dof.get_dof(x, 3), kqq[i]);
        res(dof.get_dof(x, 3)) = rq[i];
        i++;
    }

    i = 0;
    for (auto e: mesh.cathode_elements) {
        MatrixXd e_p(n, 1);
        MatrixXd e_c = MatrixXd::Zero(n, 1);
        MatrixXd e_s(n, 1);
        MatrixXd e_q(n, 1);
        MatrixXd e_t(n, 1);
        for (int j = 0; j < n; j++) {
            std::size_t node_id = mesh.elements[e * n + j];
            e_p(j, 0) = u(dof.get_dof(node_id, 0), 0);
            if (!settings::is_solid_battery) e_c(j, 0) = u(dof.get_dof(node_id, 1), 0);
            e_s(j, 0) = u(dof.get_dof(node_id, 2), 0);
            e_q(j, 0) = u(dof.get_dof(node_id, 3), 0);
            if constexpr (use_temp) {
                e_t(j, 0) = u(dof.get_dof(node_id, 4), 0);
            } else {
                e_t(j, 0) = constant::T;
            }
        }
        MatrixXd e_dc = MatrixXd::Zero(n, 1);
        MatrixXd e_ds(n, 1);
        MatrixXd e_dt(n, 1);
        for (int j = 0; j < n; j++) {
            std::size_t node_id = mesh.elements[e * n + j];
            if (!settings::is_solid_battery) e_dc(j, 0) = du(dof.get_dof(node_id, 1), 0);
            e_ds(j, 0) = du(dof.get_dof(node_id, 2), 0);
            if constexpr (use_temp) {
                e_dt(j, 0) = du(dof.get_dof(node_id, 4), 0);
            } else {
                e_dt(j, 0) = 0;
            }
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
        MatrixXd e_kts = MatrixXd::Zero(n, n);
        MatrixXd e_ktq = MatrixXd::Zero(n, n);

        MatrixXd e_rs = MatrixXd::Zero(n, 1);
        MatrixXd e_rc = MatrixXd::Zero(n, 1);
        MatrixXd e_rp = MatrixXd::Zero(n, 1);
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
            double ele_c_e = t_mat.value();
            t_mat = N.transpose() * e_t;
            double ele_c_t = t_mat.value();

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
            if (!settings::is_solid_battery) {
                double eff_2 = a * constant::l_ref * constant::l_ref * (1 - constant::trans) / d_ref;
                double eff_3 = 1 / dt * constant::l_ref * constant::l_ref / d_ref;
                e_kcc += epsilon * eff_3 * NNT * w(j) * det + d_eff * dNdNT * w(j) * det;
                e_kcq += -eff_2 * NNT * w(j) * det * j_ref / ce_int;
                if (!is_first_step) {
                    e_rc += epsilon * eff_3 * NNT * e_dc * w(j) * det + d_eff * dNdNT * e_c * w(j) * det - eff_2 * NNT * e_q *
                            w(j) * det * j_ref / ce_int;
                } else {
                    e_rc += epsilon * eff_3 * NNT * e_dc * w(j) * det + d_eff * dNdNT * e_c * w(j) * det;
                }
            }
            //e_kcc = MatrixXd::Identity(2, 2);

            // p part
            double eff_4 = a * F * constant::l_ref * constant::l_ref / k_ref;
            if (!settings::is_solid_battery) {
                e_kpp += k_eff * dNdNT * w(j) * det;
                e_kpc += d_k_eff * dNdNT * e_p * N_T * w(j) * det
                        - kd_eff / ele_c_e * dNdNT * w(j) * det
                        - d_kd_eff / ele_c_e * dNdNT * e_c * N_T * w(j) * det
                        + kd_eff / (ele_c_e * ele_c_e) * dNdNT * e_c * N_T * w(j) * det;
                e_kpq += -eff_4 * NNT * w(j) * det * j_ref;
                e_rp += k_eff * dNdNT * e_p * w(j) * det - kd_eff / ele_c_e * dNdNT * e_c * w(j) * det - eff_4 * NNT * e_q *
                        w(j) * det * j_ref;
            } else {
                e_kpp += k_eff * dNdNT * w(j) * det;
                e_kpq += -eff_4 * NNT * w(j) * det * j_ref;
                e_rp += k_eff * dNdNT * e_p * w(j) * det - eff_4 * NNT * e_q * w(j) * det * j_ref;

            }
            //e_kpp = MatrixXd::Identity(2, 2);

            // t part
            if constexpr (use_temp) {
                MatrixXd e_eta(n, 1);
                for (int l = 0; l < n; l++) {
                    e_eta(l, 0) = arr_eta[mesh.node_to_idx[mesh.elements[e * n + l]]];
                }
                double Ne_eta = (N_T * e_eta).value();
                MatrixXd e_du(n, 1);
                for (int l = 0; l < n; l++) {
                    e_du(l, 0) = arr_du[mesh.node_to_idx[mesh.elements[e * n + l]]];
                }
                // Heat transfer
                e_ktt += rho * cap * constant::l_ref * constant::l_ref * NNT / dt * w(j) * det
                        + lambda * dNdNT * w(j) * det;
                e_rt += rho * cap * constant::l_ref * constant::l_ref * NNT * e_dt / dt * w(j) * det
                        + lambda * dNdNT * e_t * w(j) * det;
                if (!is_first_step) {
                    VectorXd grad_p = dN_T * e_p;
                    VectorXd grad_s = dN_T * e_s;
                    VectorXd grad_c = dN_T * e_c;
                    double Ne_q = (N_T * e_q).value();
                    double Ne_t = (N_T * e_t).value();
                    double Ne_du = (N_T * e_du).value();

                    // Ohmic heat
                    if (!settings::is_solid_battery) e_ktt += (dk_dt * NNT / ele_c_e * grad_c.dot(grad_p)) * k_ref * w(j) * det;
                    if (!settings::is_solid_battery) e_ktp += -(k_eff * N * 2 * grad_p.transpose() * dN_T - dk_dt * N * Ne_t * grad_c.transpose() * dN_T / ele_c_e) * k_ref * w(j) * det;
                    else e_ktp += -(k_eff * N * 2 * grad_p.transpose() * dN_T) * k_ref * w(j) * det;
                    if (!settings::is_solid_battery) e_ktc += (dk_dt * N * Ne_t * grad_p.transpose() * dN_T / ele_c_e) * k_ref * w(j) * det;
                    e_kts += -(sigma_eff * N * 2 * grad_s.transpose() * dN_T) * sigma_ref * w(j) * det;
                    e_rt += -(k_eff * N * grad_p.dot(grad_p)) * k_ref * w(j) * det;
                    e_rt += -(sigma_eff * N * grad_s.dot(grad_s)) * sigma_ref * w(j) * det;
                    if (!settings::is_solid_battery) e_rt += (dk_dt * N * Ne_t / ele_c_e * grad_c.dot(grad_p)) * k_ref * w(j) * det;
                    // Irreversible heat
                    e_ktq += -(constant::F * a * NNT * Ne_eta) * constant::l_ref * constant::l_ref * j_ref * w(j) * det;
                    e_rt += -(constant::F * a * N * Ne_q * Ne_eta) * constant::l_ref * constant::l_ref * j_ref * w(j) * det;
                    // Reversible heat
                    e_ktt += -(constant::F * a * NNT * Ne_du * Ne_q) * constant::l_ref * constant::l_ref * j_ref * w(j) * det;
                    e_ktq += -(constant::F * a * NNT * Ne_du * Ne_t) * constant::l_ref * constant::l_ref * j_ref * w(j) * det;
                    e_rt += -(constant::F * a * N * Ne_du * Ne_t * Ne_q) * constant::l_ref * constant::l_ref * j_ref * w(j) * det;
                }
            }
        }

        for (int j = 0; j < n; j++) {
            std::size_t id_l = mesh.elements[e * n + j];
            for (int l = 0; l < n; l++) {
                std::size_t id_r = mesh.elements[e * n + l];
                t.emplace_back(dof.get_dof(id_l, 0), dof.get_dof(id_r, 0), e_kpp(j, l));
                if (!settings::is_solid_battery) t.emplace_back(dof.get_dof(id_l, 0), dof.get_dof(id_r, 1), e_kpc(j, l));
                t.emplace_back(dof.get_dof(id_l, 0), dof.get_dof(id_r, 3), e_kpq(j, l));

                if (!settings::is_solid_battery) {
                    t.emplace_back(dof.get_dof(id_l, 1), dof.get_dof(id_r, 1), e_kcc(j, l));
                    if (!is_first_step) t.emplace_back(dof.get_dof(id_l, 1), dof.get_dof(id_r, 3), e_kcq(j, l));
                }

                t.emplace_back(dof.get_dof(id_l, 2), dof.get_dof(id_r, 2), e_kss(j, l));
                t.emplace_back(dof.get_dof(id_l, 2), dof.get_dof(id_r, 3), e_ksq(j, l));

                if constexpr (use_temp) {
                    t.emplace_back(dof.get_dof(id_l, 4), dof.get_dof(id_r, 0), e_ktp(j, l));
                    if (!settings::is_solid_battery) t.emplace_back(dof.get_dof(id_l, 4), dof.get_dof(id_r, 1), e_ktc(j, l));
                    t.emplace_back(dof.get_dof(id_l, 4), dof.get_dof(id_r, 2), e_kts(j, l));
                    t.emplace_back(dof.get_dof(id_l, 4), dof.get_dof(id_r, 3), e_ktq(j, l));
                    t.emplace_back(dof.get_dof(id_l, 4), dof.get_dof(id_r, 4), e_ktt(j, l));
                }
            }
            res(dof.get_dof(id_l, 0)) += e_rp(j);
            if (!settings::is_solid_battery) res(dof.get_dof(id_l, 1)) += e_rc(j);
            res(dof.get_dof(id_l, 2)) += e_rs(j);
            if constexpr (use_temp) res(dof.get_dof(id_l, 4)) += e_rt(j);
        }
        i++;
    }

    delete[] kqp;
    delete[] kqc;
    delete[] kqs;
    delete[] kqq;
    delete[] rq;
    delete[] arr_eta;
    delete[] arr_du;
}

template void stiffness_cathode::generate<true>(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  const std::vector<double> &, const std::vector<double> &,
                  std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                  bool, std::vector<double> &);
template void stiffness_cathode::generate<false>(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  const std::vector<double> &, const std::vector<double> &,
                  std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                  bool, std::vector<double> &);