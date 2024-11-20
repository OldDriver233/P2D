#include "stiffness_anode_collector.h"
#include <iostream>

void stiffness_anode_collector::generate(const Eigen::Ref<MatrixXd> &u, 
                               const Eigen::Ref<MatrixXd> &du,
                               const Eigen::Ref<MatrixXd> &c_s,
                               std::vector<Eigen::Triplet<double>> &t, 
                               Eigen::Ref<VectorXd> res, 
                               bool is_first_step) {
    const double I_app = constant::I_app;
    const int dim = 1, n = 2;
    const double rho = 8940, cap = 385, lambda = 401, sigma = 5.96e6;
    int dof_cnt = this->surface_ca_coll - this->surface_an_coll + 1;
    int dof_cnt_eff = dof_cnt - (this->surface_ca_sep - this->surface_an_sep - 1);
    int dof_cnt_temp = this->points.size();
    MatrixXd xs = get_integration_point<dim, n>();
    MatrixXd w = get_integration_weight<dim, n>();

    for(int i = 0; i < this->surface_an_coll; i++) {
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

            e_ktt += rho * cap * NNT / constant::dt * w(j) * det + lambda * dNdNT * w(j) * det;
            //std::cout<<e_dt<<std::endl;
            //std::cout<<rho * cap * NNT * e_dt / constant::dt * w(j) * det<<";"<<lambda * dNdNT * e_t * w(j) * det<<std::endl;
            e_rt += rho * cap * NNT * e_dt / constant::dt * w(j) * det + lambda * dNdNT * e_t * w(j) * det - I_app * I_app / sigma * N * w(j) * det;
        }

        for(int j = 0; j < n; j++) {
            for(int l = 0; l < n; l++) {
                t.push_back(Eigen::Triplet<double>(i + j + 2 * dof_cnt + 2 * dof_cnt_eff, i + l + 2 * dof_cnt + 2 * dof_cnt_eff, e_ktt(j, l)));
            }
            res(i + j + 2 * dof_cnt + 2 * dof_cnt_eff) += e_rt(j);
        }
    }
}