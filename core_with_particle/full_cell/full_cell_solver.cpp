#include "full_cell_solver.h"
#include "../functions/functions.h"
#include <cstdlib>
#include <eigen3/Eigen/Sparse>
#include <cstdio>
#include <iostream>

inline double clamp(double x, double lower, double upper) {
    return x < lower ? lower : (x > upper ? upper : x);
}

void full_cell_solver::apply_boundary(Eigen::Ref<MatrixXd> u, Eigen::SparseMatrix<double> &k, Eigen::Ref<VectorXd> res, bool is_first_step) {
    long point_size = this->point_coord.size();
    long eff_size = point_size - (ca - an - 1);
    if(is_first_step) {
        //k.insert(0, 0) = 1;
        //k.insert(2 * point_size, 2 * point_size) = 1;
        //k.insert(2 * point_size + eff_size - 1, 2 * point_size + eff_size - 1) = 1;
        //res(0, 0) = -(0 - u(0, 0));
        //res(2 * point_size, 0) = -(uoc<1>(constant::c_int_an / constant::c_max_an) - u(2 * point_size, 0));
        //res(2 * point_size + eff_size - 1, 0) = -(uoc<2>(constant::c_int_ca / constant::c_max_ca) - u(2 * point_size + eff_size - 1, 0));
        k.insert(2 * point_size, 2 * point_size) = 1;
        res(2 * point_size, 0) = -(0 - u(2 * point_size, 0));
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

void full_cell_solver::calc(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> c_s) {
    long point_size = this->point_coord.size();
    long eff_size = point_size - (ca - an - 1);
    long element_cnt = this->point_coord.size() - 1;
    int iter_time = 0;
    double res_norm = 999999.0;
    double first_norm;
    MatrixXd du = MatrixXd::Zero(2 * point_size + 2 * eff_size, 1);
    std::vector<Eigen::Triplet<double>> coeff;

    anode_particle.pre_calc(c_s);
    cathode_particle.pre_calc(c_s);

    while(iter_time < iter && res_norm > tolerance) {
        Eigen::SparseMatrix<double> k(2 * point_size + 2 * eff_size, 2 * point_size + 2 * eff_size);
        VectorXd res = VectorXd::Zero(2 * point_size + 2 * eff_size);
        coeff.clear();
        coeff.reserve(12 * point_size);

        this->anode.generate(u, du, c_s, coeff, res, step == 0);
        this->sep.generate(u, du, c_s, coeff, res, step == 0);
        this->cathode.generate(u, du, c_s, coeff, res, step == 0);
        k.setFromTriplets(coeff.begin(), coeff.end());
        apply_boundary(u, k, res, step == 0);

        solver.compute(k);
        //std::cout<<k<<std::endl;
        MatrixXd delta = -solver.solve(res);
        //std::cout<<delta<<std::endl;
        du += delta;
        u += delta;
        anode_particle.calc(c_s, u, point_size, an, ca, 1);
        cathode_particle.calc(c_s, u, point_size, an, ca, 2);
        double norm = delta.norm();
        res_norm = res.norm() / (2 * point_size + 2 * eff_size);
        printf("Step %d Iter %d: %.12lf, %.12lf\n", step, iter_time, norm, res_norm);
        

        iter_time++;
    }
    step++;
}