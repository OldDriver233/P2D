#include "full_cell_solver.h"
#include "../functions/functions.h"
#include <cstdlib>
#include <eigen3/Eigen/Sparse>
#include <iostream>
#include <cstdio>
#include <eigen3/Eigen/SparseLU>

inline double clamp(double x, double lower, double upper) {
    return x < lower ? lower : (x > upper ? upper : x);
}

void full_cell_solver::apply_boundary(Eigen::Ref<MatrixXd> u, Eigen::SparseMatrix<double> &k, Eigen::Ref<VectorXd> res, bool is_first_step) {
    long point_size = this->point_coord.size();
    long eff_size = point_size - (ca - an - 1);
    if(is_first_step) {
        k.insert(0, 0) = 1;
        k.insert(2 * point_size, 2 * point_size) = 1;
        k.insert(2 * point_size + eff_size - 1, 2 * point_size + eff_size - 1) = 1;
        res(0, 0) = -(0 - u(0, 0));
        res(2 * point_size, 0) = -(uoc(constant::c_int_an / constant::c_max_an, 1) - u(2 * point_size, 0));
        res(2 * point_size + eff_size - 1, 0) = -(uoc(constant::c_int_ca / constant::c_max_ca, 2) - u(2 * point_size + eff_size - 1, 0));
    } else {
        k.insert(2 * point_size, 2 * point_size) = 1;
        res(2 * point_size, 0) = -(0 - u(2 * point_size, 0));

        double eff_mat_s_ca = std::pow(constant::epsilon_s_ca, constant::bruggeman);
        //double eff_mat_s_an = std::pow(constant::epsilon_s_an, constant::bruggeman);
        //double sigma_ref_an = constant::sigma_an * eff_mat_s_an;
        double sigma_ref_ca = constant::sigma_ca * eff_mat_s_ca;
        //res(2 * point_size, 0) -= 30 * constant::l_ref / sigma_ref_an;
        res(2 * point_size + eff_size - 1, 0) += 30 * constant::l_ref / sigma_ref_ca;
    }
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

        this->anode.generate(u, du, coeff, res, step == 0);
        this->sep.generate(u, du, coeff, res, step == 0);
        this->cathode.generate(u, du, coeff, res, step == 0);
        k.setFromTriplets(coeff.begin(), coeff.end());
        apply_boundary(u, k, res, step == 0);
        //std::cout<<k<<"\n"<<res<<"\n-="<<std::endl;

        solver.analyzePattern(k);
        solver.factorize(k);
        if(solver.info() != 0) {
            std::cerr<<solver.info();
        }        
        solver.compute(k);
        MatrixXd delta = - solver.solve(res);
        du += delta;
        u += delta;
        double norm = delta.norm();
        res_norm = res.norm() / (2 * point_size + 4 * eff_size);
        //std::cout<<delta<<"-"<<std::endl;
        //std::cout<<"Iter "<<iter_time<<": "<<norm<<","<<res_norm<<std::endl;
        printf("Step %d Iter %d: %lf, %lf\n", step, iter_time, norm, res_norm);
        //std::cout<<"Cond: "<<k.norm() * k.inverse().norm()<<std::endl;
        //std::cout<<k * delta + res<<"\n-"<<std::endl;
        //if(step == 42) std::cout<<u<<"-"<<std::endl;

        iter_time++;
    }
    step++;
}