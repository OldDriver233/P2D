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
}
