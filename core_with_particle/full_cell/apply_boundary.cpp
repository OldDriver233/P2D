#include "full_cell_solver.h"

void full_cell_solver::apply_boundary(Eigen::Ref<MatrixXd> u, Eigen::SparseMatrix<double> &K, Eigen::Ref<VectorXd> res,
                                      bool is_first_step) {
    const int dim = get_dim(mesh.p_type), n = get_nodes(mesh.p_type);
    for (auto x: mesh.anode_wall_nodes) {
        std::size_t dof_id = dof.get_dof(x, 2);
        K.insert(dof_id, dof_id) = 1;
        res(dof_id) = -(0 - u(dof_id, 0));
    }

    double eff_mat_s_ca = std::pow(constant::epsilon_s_ca, constant::bruggeman);
    double eff_mat_s_an = std::pow(constant::epsilon_s_an, constant::bruggeman);
    double sigma_ref_an = constant::sigma_an * eff_mat_s_an;
    double sigma_ref_ca = constant::sigma_ca * eff_mat_s_ca;

    if (dim == 1) {
        for (auto x: mesh.cathode_wall_nodes) {
            res(dof.get_dof(x, 2)) += constant::I_app * constant::l_ref / sigma_ref_ca;
        }
    } else if (dim == 2) {
        MatrixXd xs = get_integration_point<1, 2>();
        MatrixXd w = get_integration_weight<1, 2>();
        int elem_count = mesh.cathode_wall.size() / 2;

        for (int i = 0; i < elem_count; i++) {
            MatrixXd coord(2, 2);

            VectorXd i_rel = VectorXd::Zero(2);
            for (int j = 0; j < 2; j++) {
                coord.col(j) = mesh.coord.col(mesh.cathode_wall[i * 2 + j]);
            }
            for (int j = 0; j < 2; j++) {
                MatrixXd N = get_shape_func_at<1, 2>(xs.row(j).transpose());
                MatrixXd dNdu = get_shape_deriv_at<1, 2>(xs.row(j).transpose());

                MatrixXd dN = coord * dNdu;

                double ds = std::hypot(dN(0, 0), dN(1, 0));

                i_rel += constant::I_app * constant::l_ref / sigma_ref_ca * N * w(j) * ds;
            }

            for (int j = 0; j < 2; j++) {
                res(dof.get_dof(mesh.cathode_wall[i * 2 + j], 2)) += i_rel(j);
            }
        }
    }
    if (settings::calc_temperature) {
        const double t_exchange = 1;
        if (dim == 1) {
            for (auto x: mesh.anode_wall_nodes) {
                int idx = dof.get_dof(x, 4);
                K.coeffRef(idx, idx) += t_exchange * constant::l_ref;
                res(idx) += t_exchange * (u(idx, 0) - constant::t_ref) * constant::l_ref;
            }
            for (auto x: mesh.cathode_wall_nodes) {
                int idx = dof.get_dof(x, 4);
                K.coeffRef(idx, idx) += t_exchange * constant::l_ref;
                res(idx) += t_exchange * (u(idx, 0) - constant::t_ref) * constant::l_ref;
            }
        } else if (dim == 2) {
            MatrixXd xs = get_integration_point<1, 2>();
            MatrixXd w = get_integration_weight<1, 2>();

            int elem_count = mesh.cathode_wall.size() / 2;
            for (int i = 0; i < elem_count; i++) {
                MatrixXd coord(2, 2);
                VectorXd e_ref_t(2);

                MatrixXd k_ref = MatrixXd::Zero(2, 2);
                VectorXd q_rel = VectorXd::Zero(2);
                for (int j = 0; j < 2; j++) {
                    coord.col(j) = mesh.coord.col(mesh.cathode_wall[i * 2 + j]);
                    e_ref_t(j) = u(dof.get_dof(mesh.cathode_wall[i * 2 + j], 4), 0) - constant::t_ref;
                }
                for (int j = 0; j < 2; j++) {
                    MatrixXd N = get_shape_func_at<1, 2>(xs.row(j).transpose());
                    MatrixXd dNdu = get_shape_deriv_at<1, 2>(xs.row(j).transpose());

                    MatrixXd dN = coord * dNdu;

                    double ds = std::hypot(dN(0, 0), dN(1, 0));
                    double ref_t = (N.transpose() * e_ref_t).value();

                    k_ref += t_exchange * N * N.transpose() * constant::l_ref * w(j) * ds;
                    q_rel += t_exchange * N * ref_t * constant::l_ref * w(j) * ds;
                }

                for (int j = 0; j < 2; j++) {
                    int id_l = dof.get_dof(mesh.cathode_wall[i * 2 + j], 4);
                    for (int l = 0; l < 2; l++) {
                        int id_r = dof.get_dof(mesh.cathode_wall[i * 2 + l], 4);
                        K.coeffRef(id_l, id_r) += k_ref(j, l);
                    }
                    res(id_l) += q_rel(j);
                }
            }

            elem_count = mesh.anode_wall.size() / 2;
            for (int i = 0; i < elem_count; i++) {
                MatrixXd coord(2, 2);
                VectorXd e_ref_t(2);

                MatrixXd k_ref = MatrixXd::Zero(2, 2);
                VectorXd q_rel = VectorXd::Zero(2);
                for (int j = 0; j < 2; j++) {
                    coord.col(j) = mesh.coord.col(mesh.anode_wall[i * 2 + j]);
                    e_ref_t(j) = u(dof.get_dof(mesh.anode_wall[i * 2 + j], 4), 0) - constant::t_ref;
                }
                for (int j = 0; j < 2; j++) {
                    MatrixXd N = get_shape_func_at<1, 2>(xs.row(j).transpose());
                    MatrixXd dNdu = get_shape_deriv_at<1, 2>(xs.row(j).transpose());

                    MatrixXd dN = coord * dNdu;

                    double ds = std::hypot(dN(0, 0), dN(1, 0));
                    double ref_t = (N.transpose() * e_ref_t).value();

                    k_ref += t_exchange * N * N.transpose() * constant::l_ref * w(j) * ds;
                    q_rel += t_exchange * N * ref_t * constant::l_ref * w(j) * ds;
                }

                for (int j = 0; j < 2; j++) {
                    int id_l = dof.get_dof(mesh.cathode_wall[i * 2 + j], 4);
                    for (int l = 0; l < 2; l++) {
                        int id_r = dof.get_dof(mesh.cathode_wall[i * 2 + l], 4);
                        K.coeffRef(id_l, id_r) += k_ref(j, l);
                    }
                    res(id_l) += q_rel(j);
                }
            }
        }
    }
}
