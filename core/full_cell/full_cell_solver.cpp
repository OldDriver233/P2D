#include "full_cell_solver.h"
#include <eigen3/Eigen/Sparse>
#include <iostream>
#include <cstdio>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/PardisoSupport>

inline double clamp(double x, double lower, double upper) {
    return x < lower ? lower : (x > upper ? upper : x);
}

void full_cell_solver::apply_boundary(Eigen::Ref<MatrixXd> u, Eigen::SparseMatrix<double> &k, Eigen::Ref<VectorXd> res) {
    long point_size = this->point_coord.size();
    long eff_size = point_size - (ca - an - 1);

    for(int i = 0; i < 2 * point_size + 4 * eff_size; i++) {
        k.coeffRef(2 * point_size, i) = 0;
    }
    k.coeffRef(2 * point_size, 2 * point_size) = 1;
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
        Eigen::SparseMatrix<double> k(2 * point_size + 4 * eff_size, 2 * point_size + 4 * eff_size);
        VectorXd res = VectorXd::Zero(2 * point_size + 4 * eff_size);
        std::vector<Eigen::Triplet<double>> coeff;

        this->anode.generate(u, du, coeff, res);
        this->sep.generate(u, du, coeff, res);
        this->cathode.generate(u, du, coeff, res);
        k.setFromTriplets(coeff.begin(), coeff.end());
        apply_boundary(u, k, res);

        Eigen::PardisoLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(k);
        MatrixXd delta = - solver.solve(res);
        du += delta;
        u += delta;
        double norm = delta.norm();
        res_norm = res.norm() / (2 * point_size + 4 * eff_size);
        //std::cout<<delta<<"-"<<std::endl;
        //std::cout<<"Iter "<<iter_time<<": "<<norm<<","<<res_norm<<std::endl;
        printf("Iter %d: %lf, %lf\n", iter_time, norm, res_norm);
        //std::cout<<k<<"\n"<<res<<"\n-="<<std::endl;
        //std::cout<<"Cond: "<<k.norm() * k.inverse().norm()<<std::endl;
        //std::cout<<k * delta + res<<"\n-"<<std::endl;

        iter_time++;
    }
    step++;
}