#include "full_cell_solver.h"
#include <iostream>
#include <eigen3/Eigen/IterativeLinearSolvers>

inline double clamp(double x, double lower, double upper) {
    return x < lower ? lower : (x > upper ? upper : x);
}

void full_cell_solver::apply_boundary(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> k, Eigen::Ref<VectorXd> res) {
    long point_size = this->point_coord.size();
    long eff_size = point_size - (ca - an - 1);

    k.row(2 * point_size) = MatrixXd::Zero(1, 2 * point_size + 4 * eff_size);
    k(2 * point_size, 2 * point_size) = 1;
    res(2 * point_size, 0) = -(0 - u(0, 0));

    double eff_mat_s = std::pow(constant::epsilon_s_ca, constant::bruggeman);
    double sigma_ref = constant::sigma_ca * eff_mat_s;
    res(2 * point_size + eff_size - 1, 0) += 30 * constant::l_ref / sigma_ref;
}

void full_cell_solver::calc(Eigen::Ref<MatrixXd> u) {
    long point_size = this->point_coord.size();
    long eff_size = point_size - (ca - an - 1);
    long element_cnt = this->point_coord.size() - 1;
    int iter_time = 0;
    double res_norm = 999999.0;
    double first_norm;
    MatrixXd du = MatrixXd::Zero(2 * point_size + 4 * eff_size, 1);

    while(iter_time < iter && res_norm > tolerance) {
        MatrixXd k = MatrixXd::Zero(2 * point_size + 4 * eff_size, 2 * point_size + 4 * eff_size);
        VectorXd res = VectorXd::Zero(2 * point_size + 4 * eff_size);

        this->anode.generate(u, du, k, res);
        this->sep.generate(u, du, k, res);
        this->cathode.generate(u, du, k, res);
        apply_boundary(u, k, res);

        MatrixXd delta = - k.lu().solve(res);
        du += delta;
        u += delta;
        double norm = delta.norm();
        res_norm = res.norm() / (2 * point_size + 4 * eff_size);
        //std::cout<<delta<<"-"<<std::endl;
        std::cout<<"Iter "<<iter_time<<": "<<norm<<","<<res_norm<<std::endl;
        //std::cout<<k<<"\n"<<res<<"\n-="<<std::endl;
        //std::cout<<"Cond: "<<k.norm() * k.inverse().norm()<<std::endl;
        //std::cout<<k * delta + res<<"\n-"<<std::endl;

        iter_time++;
    }
    step++;
}