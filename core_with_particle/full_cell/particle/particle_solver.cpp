#include "particle_solver.h"

#include <ranges>
#include "../../io/settings/settings.h"
#include "../../step_control/StepControl.h"

void particle_solver::pre_calc(const Eigen::Ref<MatrixXd> &c_s) {
    this->last_cs = c_s;
}

double particle_solver::calc(Eigen::Ref<MatrixXd> c_s, const Eigen::Ref<MatrixXd> &u, int type, const mesh_reader& mesh, const dof_assigner& dof) {
    int particle_dof_size = constant::particle_segment + 1;
    double ret = 0.0;
    if (type == 1) {
        int i = 0;
        for (auto x: dof.particle_to_node) {
            if (mesh.anode_nodes.contains(x)) {
                if (settings::calc_temperature || settings::use_adaptive_time_step) {
                    double T;
                    if (settings::calc_temperature) {
                        T = u(dof.get_dof(x, 4), 0);
                    } else {
                        T = constant::t_ref;
                    }
                    double coeff = exp(constant::anode.diffuse_energy / constant::R * (-1 / T + 1 / constant::t_ref));
                    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_A;
                    solver_A.compute(assembled_A_1 * coeff + assembled_A_2 / this->step_control->dt_now);
                    j_coeff = solver_A.solve(pre_j_coeff);
                    MatrixXd Bc_s = assembled_B / this->step_control->dt_now * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
                    c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) =
                    solver_A.solve(Bc_s) - j_coeff * u(dof.get_dof(x, 3), 0);
                    ret = j_coeff(constant::particle_segment);
                } else {
                    MatrixXd Bc_s = assembled_B / this->step_control->dt_now * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
                    c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) =
                    constant_solver.solve(Bc_s) - j_coeff * u(dof.get_dof(x, 3), 0);
                }
            }
            i++;
        }
    } else {
        int i = 0;
        for (auto x: dof.particle_to_node) {
            if (mesh.cathode_nodes.contains(x)) {
                if (settings::calc_temperature || settings::use_adaptive_time_step) {
                    double T;
                    if (settings::calc_temperature) {
                        T = u(dof.get_dof(x, 4), 0);
                    } else {
                        T = constant::t_ref;
                    }
                    double coeff = exp(constant::cathode.diffuse_energy / constant::R * (-1 / T + 1 / constant::t_ref));
                    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_A;
                    solver_A.compute(assembled_A_1 * coeff + assembled_A_2 / this->step_control->dt_now);
                    j_coeff = solver_A.solve(pre_j_coeff);
                    MatrixXd Bc_s = assembled_B / this->step_control->dt_now * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
                    c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) =
                    solver_A.solve(Bc_s) - j_coeff * u(dof.get_dof(x, 3), 0);
                    ret = j_coeff(constant::particle_segment);
                } else {
                    MatrixXd Bc_s = assembled_B / this->step_control->dt_now * last_cs.block(i * particle_dof_size, 0, particle_dof_size, 1);
                    c_s.block(i * particle_dof_size, 0, particle_dof_size, 1) =
                    constant_solver.solve(Bc_s) - j_coeff * u(dof.get_dof(x, 3), 0);
                }
            }
            i++;
        }
    }
    return ret;
}

double particle_solver::calc_stress(Eigen::Ref<MatrixXd> c_s, const Eigen::Ref<MatrixXd> &u, int type, const mesh_reader &mesh, const dof_assigner &dof) {
    int particle_dof_size = constant::particle_segment + 1;
    double ret = 0.0;
    const int dim = 1, n = 2;
    double dt;
    if (settings::use_adaptive_time_step) {
        dt = step_control->dt_now;
    } else {
        dt = constant::dt;
    }
    MatrixXd xs = get_integration_point<dim, n>();
    MatrixXd w = get_integration_weight<dim, n>();

    if (type == 1) {
        int idx = 0;
        for (auto x: dof.particle_to_node) {
            std::vector<Eigen::Triplet<double>> k_coeff;
            Eigen::SparseMatrix<double> K(particle_dof_size, particle_dof_size);
            VectorXd b;

            if (mesh.anode_nodes.contains(x)) {
                double T;
                double omega = constant::anode.omega, E = constant::anode.E, nu = constant::anode.nu;
                if (settings::calc_temperature) {
                    T = u(dof.get_dof(x, 4), 0);
                } else {
                    T = constant::t_ref;
                }
                double Z = 2 * omega * omega * E / (9 * constant::R * T * (1 - nu));
                double D = D_s * exp(constant::anode.diffuse_energy / constant::R * (-1 / T + 1 / constant::t_ref));

                int iter_round = 0;
                double norm = 1, rel_norm = 1;
                while (iter_round < 1 && norm / rel_norm > 1e-3) {
                    K.setZero();
                    b = VectorXd::Zero(particle_dof_size);
                    k_coeff.clear();

                    for (int i = 0; i < constant::particle_segment; i++) {
                        MatrixXd coords(1, n);
                        coords << point_coord(i), point_coord(i + 1);
                        MatrixXd e_ss(n, 1);
                        e_ss << c_s(idx * particle_dof_size + i, 0), c_s(idx * particle_dof_size + i + 1, 0);
                        MatrixXd e_dss(n, 1);
                        e_dss << last_cs(idx * particle_dof_size + i, 0), last_cs(idx * particle_dof_size + i + 1, 0);
                        e_dss = e_ss - e_dss;

                        MatrixXd e_k = MatrixXd::Zero(n, n);
                        VectorXd e_b = VectorXd::Zero(n);

                        for (int j = 0; j < n; j++) {
                            MatrixXd N = cached_N[i * n + j];
                            double det_J = cached_J[i * n + j];
                            MatrixXd dN = cached_dN[i * n + j];
                            MatrixXd N_T = N.transpose();
                            MatrixXd dN_T = dN.transpose();

                            double s = xs(j);
                            double lower = coords(0, 0);
                            double upper = coords(0, 1);
                            double x = lower + (s + 1) * (upper - lower) / 2;
                            double ele_c = (N_T * e_ss).value();

                            double eff_1 = R_s * R_s / D_sref;
                            double eff_2 = D / D_sref;
                            double eff_3 = D / D_sref * Z * c_max;

                            e_k += N * N_T * x * x * w(j) * det_J * eff_1 / dt;
                            e_b += N * N_T * x * x * w(j) * det_J * eff_1 * e_dss / dt;
                            e_k += dN * dN_T * x * x * w(j) * det_J * eff_2;
                            e_b += dN * dN_T * x * x * w(j) * det_J * eff_2 * e_ss;
                            e_k += dN * dN_T * x * x * ele_c * w(j) * det_J * eff_3;
                            e_b += dN * dN_T * x * x * ele_c * w(j) * det_J * eff_3 * e_ss;
                        }

                        for (int j = 0; j < n; j++) {
                            for (int l = 0; l < n; l++) {
                                k_coeff.emplace_back(i + j, i + l, e_k(j, l));
                            }
                            b(i + j) += e_b(j);
                        }
                    }
                    K.setFromTriplets(k_coeff.begin(), k_coeff.end());
                    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
                    solver.compute(K);
                    j_coeff = solver.solve(pre_j_coeff);
                    ret = j_coeff(constant::particle_segment);

                    b(particle_dof_size - 1) += R_s / D_sref * constant::j_ref / c_max * u(dof.get_dof(x, 3), 0);
                    MatrixXd delta = -solver.solve(b);
                    c_s.block(idx * particle_dof_size, 0, particle_dof_size, 1) += delta;
                    if (iter_round == 0) rel_norm = delta.norm();
                    norm = delta.norm();
                    iter_round++;

                    //std::cout<<norm<<std::endl;
                    //std::cout<<norm / rel_norm<<std::endl;
                    //std::cout<<"\n";
                }
            }
            idx++;
        }
    } else {
        int idx = 0;
        for (auto x: dof.particle_to_node) {
            std::vector<Eigen::Triplet<double>> k_coeff;
            Eigen::SparseMatrix<double> K(particle_dof_size, particle_dof_size);
            VectorXd b;

            if (mesh.cathode_nodes.contains(x)) {
                double T;
                double omega = constant::cathode.omega, E = constant::cathode.E, nu = constant::cathode.nu;
                if (settings::calc_temperature) {
                    T = u(dof.get_dof(x, 4), 0);
                } else {
                    T = constant::t_ref;
                }
                double Z = 2 * omega * omega * E / (9 * constant::R * T * (1 - nu));
                double D = D_s * exp(constant::cathode.diffuse_energy / constant::R * (-1 / T + 1 / constant::t_ref));

                int iter_round = 0;
                double norm = 1, rel_norm = 1;
                while (iter_round < 1 && norm / rel_norm > 1e-3) {
                    K.setZero();
                    b = VectorXd::Zero(particle_dof_size);
                    k_coeff.clear();

                    for (int i = 0; i < constant::particle_segment; i++) {
                        MatrixXd coords(1, n);
                        coords << point_coord(i), point_coord(i + 1);
                        MatrixXd e_ss(n, 1);
                        e_ss << c_s(idx * particle_dof_size + i, 0), c_s(idx * particle_dof_size + i + 1, 0);
                        MatrixXd e_dss(n, 1);
                        e_dss << last_cs(idx * particle_dof_size + i, 0), last_cs(idx * particle_dof_size + i + 1, 0);
                        e_dss = e_ss - e_dss;

                        MatrixXd e_k = MatrixXd::Zero(n, n);
                        VectorXd e_b = VectorXd::Zero(n);

                        for (int j = 0; j < n; j++) {
                            MatrixXd N = cached_N[i * n + j];
                            double det_J = cached_J[i * n + j];
                            MatrixXd dN = cached_dN[i * n + j];
                            MatrixXd N_T = N.transpose();
                            MatrixXd dN_T = dN.transpose();

                            double s = xs(j);
                            double lower = coords(0, 0);
                            double upper = coords(0, 1);
                            double x = lower + (s + 1) * (upper - lower) / 2;
                            double ele_c = (N_T * e_ss).value();

                            double eff_1 = R_s * R_s / D_sref;
                            double eff_2 = D / D_sref;
                            double eff_3 = D / D_sref * Z * c_max;

                            e_k += N * N_T * x * x * w(j) * det_J * eff_1 / dt;
                            e_b += N * N_T * x * x * w(j) * det_J * eff_1 * e_dss / dt;
                            e_k += dN * dN_T * x * x * w(j) * det_J * eff_2;
                            e_b += dN * dN_T * x * x * w(j) * det_J * eff_2 * e_ss;
                            e_k += dN * dN_T * x * x * ele_c * w(j) * det_J * eff_3;
                            e_b += dN * dN_T * x * x * ele_c * w(j) * det_J * eff_3 * e_ss;
                        }

                        for (int j = 0; j < n; j++) {
                            for (int l = 0; l < n; l++) {
                                k_coeff.emplace_back(i + j, i + l, e_k(j, l));
                            }
                            b(i + j) += e_b(j);
                        }
                    }
                    K.setFromTriplets(k_coeff.begin(), k_coeff.end());
                    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
                    solver.compute(K);
                    j_coeff = solver.solve(pre_j_coeff);
                    ret = j_coeff(constant::particle_segment);

                    b(particle_dof_size - 1) += R_s / D_sref * constant::j_ref / c_max * u(dof.get_dof(x, 3), 0);
                    MatrixXd delta = -solver.solve(b);
                    c_s.block(idx * particle_dof_size, 0, particle_dof_size, 1) += delta;
                    if (iter_round == 0) rel_norm = delta.norm();
                    norm = delta.norm();
                    iter_round++;
                }
            }
            idx++;
        }
    }

    return ret;
}

void particle_solver::get_average_concentration(const Eigen::Ref<MatrixXd> &c_s, std::vector<double> &conc, int type, const mesh_reader& mesh, const dof_assigner& dof) {
    const int dim = 1, n = 2;
    MatrixXd xs = get_integration_point<dim, n>();
    MatrixXd w = get_integration_weight<dim, n>();
    int particle_dof_size = constant::particle_segment + 1;
    if (type == 1) {
        int idx = 0;
        for (auto x: dof.particle_to_node) {
            if (mesh.anode_nodes.contains(x)) {
                double val = 0;
                for (int i = 0; i < constant::particle_segment; i++) {
                    MatrixXd coords(1, n);
                    coords << point_coord(i), point_coord(i + 1);
                    MatrixXd e_ss(n, 1);
                    e_ss << c_s(idx * particle_dof_size + i, 0), c_s(idx * particle_dof_size + i + 1, 0);
                    for (int j = 0; j < n; j++) {
                        MatrixXd N = get_shape_func_at<dim, n>(xs.row(j).transpose());
                        MatrixXd N_T = N.transpose();

                        double s = xs(j);
                        double lower = coords(0, 0);
                        double upper = coords(0, 1);
                        double x = lower + (s + 1) * (upper - lower) / 2;
                        double ele_c = (N_T * e_ss).value();
                        double length = (upper - lower) / 2;

                        val += x * x * ele_c * length * w(j);
                    }
                }
                conc[idx] = val * 3 * c_max;
            }
            idx++;
        }
    } else {
        int idx = 0;
        for (auto x: dof.particle_to_node) {
            if (mesh.cathode_nodes.contains(x)) {
                double val = 0;
                for (int i = 0; i < constant::particle_segment; i++) {
                    MatrixXd coords(1, n);
                    coords << point_coord(i), point_coord(i + 1);
                    MatrixXd e_ss(n, 1);
                    e_ss << c_s(idx * particle_dof_size + i, 0), c_s(idx * particle_dof_size + i + 1, 0);
                    for (int j = 0; j < n; j++) {
                        MatrixXd N = get_shape_func_at<dim, n>(xs.row(j).transpose());
                        MatrixXd N_T = N.transpose();

                        double s = xs(j);
                        double lower = coords(0, 0);
                        double upper = coords(0, 1);
                        double x = lower + (s + 1) * (upper - lower) / 2;
                        double ele_c = (N_T * e_ss).value();
                        double length = (upper - lower) / 2;

                        val += x * x * ele_c * length * w(j);
                    }
                }
                conc[idx] = val * 3 * c_max;
            }
            idx++;
        }
    }
}
