#include "full_cell_solver.h"
#include <iostream>

inline double clamp(double x, double lower, double upper) {
    return x < lower ? lower : (x > upper ? upper : x);
}

void full_cell_solver::apply_boundary(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> k, Eigen::Ref<VectorXd> res) {
    long point_size = this->point_coord.size();
    long eff_size = point_size - (ca - an - 1);

    k.row(0) = MatrixXd::Zero(1, 2 * point_size + 4 * eff_size);
    k(0, 0) = 1;
    res(0, 0) = -(0 - u(0, 0));

    //res(2 * point_size, 0) += 1e-2;
    res(2 * point_size + eff_size - 1, 0) -= 1e-2;
}

void full_cell_solver::calc(Eigen::Ref<MatrixXd> u) {
    long point_size = this->point_coord.size();
    long eff_size = point_size - (ca - an - 1);
    long element_cnt = this->point_coord.size() - 1;
    int iter_time = 0;
    double ratio = 999999.0;
    double first_norm;
    MatrixXd du = MatrixXd::Zero(2 * point_size + 4 * eff_size, 1);

    while(iter_time < iter && ratio > tolerance) {
        MatrixXd k = MatrixXd::Zero(2 * point_size + 4 * eff_size, 2 * point_size + 4 * eff_size);
        VectorXd res = VectorXd::Zero(2 * point_size + 4 * eff_size);

        this->anode.generate(u, du, k, res);
        this->sep.generate(u, du, k, res);
        this->cathode.generate(u, du, k, res);
        this->apply_boundary(u, k, res);
        MatrixXd delta = - k.lu().solve(res);
        du += delta;
        u += delta;
        for(int i = 2 * point_size + 2 * eff_size; i < 2 * point_size + 2 * eff_size + an + 1; i++) {
            u(i, 0) = clamp(u(i, 0), 0, 30555);
        }
        for(int i = 2 * point_size + 2 * eff_size + an + 1; i < 2 * point_size + 3 * eff_size; i++) {
            u(i, 0) = clamp(u(i, 0), 0, 51554);
        }
        for(int i = 2 * point_size + 3 * eff_size; i < 2 * point_size + 3 * eff_size + an + 1; i++) {
            u(i, 0) = clamp(u(i, 0), 0, 30555);
        }
        for(int i = 2 * point_size + 3 * eff_size + an + 1; i < 2 * point_size + 4 * eff_size; i++) {
            u(i, 0) = clamp(u(i, 0), 0, 51554);
        }
        double norm = delta.norm();
        if(iter_time != 0) {
            ratio = norm / first_norm;
        } else {
            first_norm = norm;
        }
        std::cout<<"Iter "<<iter_time<<": "<<ratio<<std::endl;

        iter_time++;
    }
    step++;
}