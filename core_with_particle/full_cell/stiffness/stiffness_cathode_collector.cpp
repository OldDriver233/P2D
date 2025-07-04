#include "stiffness_cathode_collector.h"

void stiffness_cathode_collector::generate(const Eigen::Ref<MatrixXd> &u, 
                                           const Eigen::Ref<MatrixXd> &du,
                                           const Eigen::Ref<MatrixXd> &c_s,
                                           std::vector<Eigen::Triplet<double>> &t,
                                           Eigen::Ref<VectorXd> res,
                                           bool is_first_step) {
    const double I_app = constant::I_app;
    const int dim = get_dim(mesh.p_type), n = get_nodes(mesh.p_type);
    std::size_t node_cnt = mesh.node_count;
    const double rho = constant::density_ca_collector, cap = constant::capacity_ca_collector, lambda = constant::lambda_ca_collector, sigma = constant::sigma_ca_collector;

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

    int i = 0;
    for (auto e: mesh.cathode_cc_elements) {
        MatrixXd e_ktt = MatrixXd::Zero(n, n);
        MatrixXd e_rt = MatrixXd::Zero(n, 1);
        MatrixXd e_t(n, 1);
        MatrixXd e_dt(n, 1);

        for (int j = 0; j < n; j++) {
            e_t(j, 0) = u(dof.get_dof(mesh.elements[e * n + j], 4), 0);
            e_dt(j, 0) = du(dof.get_dof(mesh.elements[e * n + j], 4), 0);
        }

        for (int j = 0; j < n; j++) {
            const MatrixXd &N = shapes.cached_matrix_N[e * n + j];
            const MatrixXd &dN = shapes.cached_matrix_dN[e * n + j];
            const MatrixXd &NNT = shapes.cached_matrix_NNT[e * n + j];
            const MatrixXd &dNdNT = shapes.cached_matrix_dNdNT[e * n + j];
            MatrixXd N_T = N.transpose();
            MatrixXd dN_T = dN.transpose();
            double det = shapes.cached_det_J[e * n + j];
            double s = xs(j);

            e_ktt += rho * cap * constant::l_ref * constant::l_ref * NNT / st->dt_now * w(j) * det
                     + lambda * dNdNT * w(j) * det;
            e_rt += rho * cap * constant::l_ref * constant::l_ref * NNT * e_dt / st->dt_now * w(j) * det
                    + lambda * dNdNT * e_t * w(j) * det
                    - I_app * I_app * constant::l_ref * constant::l_ref / sigma * N * w(j) * det;
        }

        for (int j = 0; j < n; j++) {
            std::size_t id_l = mesh.elements[e * n + j];
            for (int l = 0; l < n; l++) {
                std::size_t id_r = mesh.elements[e * n + l];
                t.emplace_back(dof.get_dof(id_l, 4), dof.get_dof(id_r, 4), e_ktt(j, l));
            }
            res(dof.get_dof(id_l, 4)) += e_rt(j);
        }
        i++;
    }
}
