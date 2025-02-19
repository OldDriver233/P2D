#include "full_cell_solver.h"
#include "../functions/functions.h"
#include <cstdlib>
#include <eigen3/Eigen/Sparse>
#include <cstdio>
#include <iostream>

inline double clamp(double x, double lower, double upper) {
    return x < lower ? lower : (x > upper ? upper : x);
}

void full_cell_solver::apply_boundary(Eigen::Ref<MatrixXd> u, Eigen::SparseMatrix<double> &k, Eigen::Ref<VectorXd> res,
                                      bool is_first_step) {
    long point_size = this->cacoll - this->ancoll + 1;
    long eff_size = point_size - (ca - an - 1);
    long all_size = this->point_coord.size();
    if (is_first_step) {
        //k.insert(0, 0) = 1;
        //k.insert(2 * point_size, 2 * point_size) = 1;
        //k.insert(2 * point_size + eff_size - 1, 2 * point_size + eff_size - 1) = 1;
        //res(0, 0) = -(0 - u(0, 0));
        //res(2 * point_size, 0) = -(uoc<1>(constant::c_int_an / constant::c_max_an) - u(2 * point_size, 0));
        //res(2 * point_size + eff_size - 1, 0) = -(uoc<2>(constant::c_int_ca / constant::c_max_ca) - u(2 * point_size + eff_size - 1, 0));
        k.insert(2 * point_size, 2 * point_size) = 1;
        res(2 * point_size, 0) = -(0 - u(2 * point_size, 0));
        double eff_mat_s_ca = std::pow(constant::epsilon_s_ca, constant::bruggeman);
        double eff_mat_s_an = std::pow(constant::epsilon_s_an, constant::bruggeman);
        double sigma_ref_an = constant::sigma_an * eff_mat_s_an;
        double sigma_ref_ca = constant::sigma_ca * eff_mat_s_ca;
        //res(2 * point_size + (an - ancoll), 0) -= constant::I_app * constant::l_ref / sigma_ref_an;
        //res(2 * point_size + (an - ancoll) + 1, 0) += constant::I_app * constant::l_ref / sigma_ref_ca;
        res(2 * point_size, 0) -= constant::I_app * constant::l_ref / sigma_ref_an;
        res(2 * point_size + eff_size - 1, 0) += constant::I_app * constant::l_ref / sigma_ref_ca;
        //k.insert(0, 0) = 1;
        //res(0, 0) -= (0 - u(0, 0));
    } else {
        k.insert(2 * point_size, 2 * point_size) = 1;
        res(2 * point_size, 0) = -(0 - u(2 * point_size, 0));
        //k.insert(0, 0) = 1;
        //res(0, 0) -= (0 - u(0, 0));

        double eff_mat_s_ca = std::pow(constant::epsilon_s_ca, constant::bruggeman);
        double eff_mat_s_an = std::pow(constant::epsilon_s_an, constant::bruggeman);
        double sigma_ref_an = constant::sigma_an * eff_mat_s_an;
        double sigma_ref_ca = constant::sigma_ca * eff_mat_s_ca;
        //res(2 * point_size + (an - ancoll), 0) -= constant::I_app * constant::l_ref / sigma_ref_an;
        //res(2 * point_size + (an - ancoll) + 1, 0) += constant::I_app * constant::l_ref / sigma_ref_ca;
        res(2 * point_size, 0) -= constant::I_app * constant::l_ref / sigma_ref_an;
        res(2 * point_size + eff_size - 1, 0) += constant::I_app * constant::l_ref / sigma_ref_ca;
    }
    if (settings::calc_temperature) {
        const double t_exchange = 1;
        k.coeffRef(2 * point_size + 2 * eff_size, 2 * point_size + 2 * eff_size) += t_exchange * constant::l_ref;
        k.coeffRef(2 * point_size + 2 * eff_size + all_size - 1, 2 * point_size + 2 * eff_size + all_size - 1) += t_exchange * constant::l_ref;
        res(2 * point_size + 2 * eff_size) -= t_exchange * (constant::t_ref - u(2 * point_size + 2 * eff_size, 0)) * constant::l_ref;
        res(2 * point_size + 2 * eff_size + all_size - 1) += t_exchange * (
            u(2 * point_size + 2 * eff_size + all_size - 1, 0) - constant::t_ref) * constant::l_ref;
    }
}

void full_cell_solver::calc(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> c_s) {
    long point_size = this->cacoll - this->ancoll + 1;
    long eff_size = point_size - (ca - an - 1);
    long element_cnt = point_size - 1;
    long all_size = this->point_coord.size();
    int iter_time = 0;
    double first_norm, first_delta_norm;
    double res_norm = 1.0;
    double rel_tol = 1.0;
    double rel_delta = 1.0;
    MatrixXd du;
    if (settings::calc_temperature) du = MatrixXd::Zero(2 * point_size + 2 * eff_size + all_size, 1);
    else du = MatrixXd::Zero(2 * point_size + 2 * eff_size, 1);
    std::vector<Eigen::Triplet<double> > coeff;
    std::vector<double> temp;

    anode_particle.pre_calc(c_s);
    cathode_particle.pre_calc(c_s);
    //printf("Step\tIter\tRelTol\tDelta\n");

    while (iter_time < iter && rel_tol > tolerance) {
        Eigen::SparseMatrix<double> k(2 * point_size + 2 * eff_size + all_size,
                                      2 * point_size + 2 * eff_size + all_size);
        if (!settings::calc_temperature) k.resize(2 * point_size + 2 * eff_size, 2 * point_size + 2 * eff_size);
        VectorXd res;
        if (settings::calc_temperature) res = VectorXd::Zero(2 * point_size + 2 * eff_size + all_size);
        else res = VectorXd::Zero(2 * point_size + 2 * eff_size);
        coeff.clear();
        coeff.reserve(16 * 5 * all_size);
        temp.clear();
        temp.reserve(3 * point_size);

        if (settings::calc_temperature) {
            this->anode.generate<true>(u, du, c_s, coeff, res, step == 0, temp);
            this->sep.generate<true>(u, du, c_s, coeff, res, step == 0, temp);
            this->cathode.generate<true>(u, du, c_s, coeff, res, step == 0, temp);
            this->anode_collector.generate(u, du, c_s, coeff, res, step == 0);
            this->cathode_collector.generate(u, du, c_s, coeff, res, step == 0);
        } else {
            this->anode.generate<false>(u, du, c_s, coeff, res, step == 0, temp);
            this->sep.generate<false>(u, du, c_s, coeff, res, step == 0, temp);
            this->cathode.generate<false>(u, du, c_s, coeff, res, step == 0, temp);
        }
        k.setFromTriplets(coeff.begin(), coeff.end());
        apply_boundary(u, k, res, false);

        solver.compute(k);
        //std::cout<<k<<std::endl;
        //std::cout<<res<<std::endl;
        MatrixXd delta = -solver.solve(res);
        du += delta;
        u += delta;

        if (step != 0) {
            if (settings::calc_temperature) {
                double j1, j2;
                j1 = -anode_particle.calc<false>(c_s, u, point_size, an - ancoll, ca - ancoll, 1, 2 * point_size + 2 * eff_size + ancoll);
                j2 = -cathode_particle.calc<false>(c_s, u, point_size, an - ancoll, ca - ancoll, 2, 2 * point_size + 2 * eff_size + ancoll);
                anode.dc_ssdj = j1;
                cathode.dc_ssdj = j2;
            } else {
                anode_particle.calc<true>(c_s, u, point_size, an - ancoll, ca - ancoll, 1, 2 * point_size + 2 * eff_size + ancoll);
                cathode_particle.calc<true>(c_s, u, point_size, an - ancoll, ca - ancoll, 2, 2 * point_size + 2 * eff_size + ancoll);
            }
        }
        double norm = delta.norm();
        res_norm = res.norm();
        if (iter_time == 0) {
            first_norm = res_norm;
            first_delta_norm = norm;
        } else {
            rel_tol = res_norm / first_norm;
            rel_delta = norm / first_delta_norm;
        }
        //printf("Step %d Iter %d: %.12lf, %.12lf\n", step, iter_time, norm, res_norm);
        printf("%-8d%-8d%1.5lf %1.6lf\n", step, iter_time, rel_tol, rel_delta);
        /*
        if (step == 149) {
            printf("%-8d%-8d%1.5lf %1.6lf\n", step, iter_time, rel_tol, rel_delta);
            std::cout<<res<<std::endl;
            std::cout<<std::endl;
        }
        */

        iter_time++;
    }
    if (step % 36 == 0) {
        VectorXd q_ohm = VectorXd::Zero(element_coord.size());
        VectorXd q_rxn = VectorXd::Zero(element_coord.size());
        VectorXd q_rev = VectorXd::Zero(element_coord.size());
        VectorXd v_eta = VectorXd::Zero(element_coord.size());
        for (int i = 0; i < temp.size(); i++) {
            if (i % 4 == 0) q_ohm(i / 4) = temp[i];
            if (i % 4 == 1) q_rxn(i / 4) = temp[i];
            if (i % 4 == 2) q_rev(i / 4) = temp[i];
            if (i % 4 == 3) v_eta(i / 4) = temp[i];
        }
        Q_ohm.append(q_ohm, step);
        Q_rxn.append(q_rxn, step);
        Q_rev.append(q_rev, step);
        eta.append(v_eta, step);
    }
    step++;
}
