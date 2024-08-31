#include "stiffness_anode.h"
#include "../../constants/constant.h"
#include <iostream>

void stiffness_anode::generate(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> du, Eigen::Ref<MatrixXd> k, Eigen::Ref<VectorXd> res) {
    int dim = 1, n = 2;
    int dof_cnt = this->points.size();
    int dof_cnt_eff = dof_cnt - (this->surface_ca_sep - this->surface_an_sep - 1);
    int elem_cnt = this->points.size() - 1;
    MatrixXd xs = get_integration_point(dim, n);
    MatrixXd w = get_integration_weight(dim, n);
    double dt = constant::dt;

    for(int i = 0; i < this->surface_an_sep; i++) {
        MatrixXd e_p = u({i, i + 1}, 0);
        MatrixXd e_c = u({dof_cnt + i, dof_cnt + i + 1}, 0);
    }
}