#include "full_cell_solver.h"
#include <iostream>

void full_cell_solver::calc(Eigen::Ref<MatrixXd> u) {
    long element_cnt = this->point_coord.size() - 1;
    int iter_time = 0;
    double ratio = 999999.0;
    double first_norm;
    MatrixXd du = MatrixXd::Zero(2 * element_cnt + 2, 1);

    while(iter_time < iter && ratio > tolerance) {
        MatrixXd k = MatrixXd::Zero(2 * element_cnt + 2, 2 * element_cnt + 2);
        VectorXd res = VectorXd::Zero(2 * element_cnt + 2);

        this->sep.generate(u, du, k, res);
        MatrixXd delta = - k.lu().solve(res);
        du += delta;
        u += delta;
        double norm = delta.norm();
        if(iter_time != 0) {
            ratio = norm / first_norm;
        } else {
            first_norm = norm;
        }

        iter_time++;
    }
}