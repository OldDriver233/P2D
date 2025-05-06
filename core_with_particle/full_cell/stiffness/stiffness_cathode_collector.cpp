#include "stiffness_cathode_collector.h"

void stiffness_cathode_collector::generate(const Eigen::Ref<MatrixXd> &u, 
                                           const Eigen::Ref<MatrixXd> &du,
                                           const Eigen::Ref<MatrixXd> &c_s,
                                           std::vector<Eigen::Triplet<double>> &t,
                                           Eigen::Ref<VectorXd> res,
                                           bool is_first_step) {
    const double I_app = constant::I_app;
    const int dim = 1, n = 2;
    const double rho = constant::density_ca_collector, cap = constant::capacity_ca_collector, lambda = constant::lambda_ca_collector, sigma = constant::sigma_ca_collector;
    int dof_cnt = this->surface_ca_coll - this->surface_an_coll + 1;
    int dof_cnt_eff = dof_cnt - (this->surface_ca_sep - this->surface_an_sep - 1);
    int dof_cnt_temp = this->points.size();
    MatrixXd xs = get_integration_point<dim, n>();
    MatrixXd w = get_integration_weight<dim, n>();

    for(int i = this->surface_ca_coll; i < dof_cnt_temp - 1; i++) {
        MatrixXd e_ktt = MatrixXd::Zero(n, n);
        MatrixXd e_rt = MatrixXd::Zero(n, 1);

        MatrixXd e_t = u({2 * dof_cnt + 2 * dof_cnt_eff + i, 2 * dof_cnt + 2 * dof_cnt_eff + i + 1}, 0);
        MatrixXd e_dt = du({2 * dof_cnt + 2 * dof_cnt_eff + i, 2 * dof_cnt + 2 * dof_cnt_eff + i + 1}, 0);

        for(int j = 0; j < n; j++) {
            const MatrixXd &N = cached_matrix_N[i * n + j];
            const MatrixXd &dN = cached_matrix_dN[i * n + j];
            const MatrixXd &NNT = cached_matrix_NNT[i * n + j];
            const MatrixXd &dNdNT = cached_matrix_dNdNT[i * n + j];
            MatrixXd N_T = N.transpose();
            MatrixXd dN_T = dN.transpose();
            double det = cached_det_J[i * n + j];
            double s = xs(j);

            e_ktt += rho * cap * constant::l_ref * constant::l_ref * NNT / st->dt_now * w(j) * det
                     + lambda * dNdNT * w(j) * det;
            e_rt += rho * cap * constant::l_ref * constant::l_ref * NNT * e_dt / st->dt_now * w(j) * det
                    + lambda * dNdNT * e_t * w(j) * det 
                    - I_app * I_app * constant::l_ref * constant::l_ref / sigma * N * w(j) * det;
        }

        for(int j = 0; j < n; j++) {
            for(int l = 0; l < n; l++) {
                t.emplace_back(i + j + 2 * dof_cnt + 2 * dof_cnt_eff, i + l + 2 * dof_cnt + 2 * dof_cnt_eff, e_ktt(j, l));
            }
            res(i + j + 2 * dof_cnt + 2 * dof_cnt_eff) += e_rt(j);
        }
    }
}
