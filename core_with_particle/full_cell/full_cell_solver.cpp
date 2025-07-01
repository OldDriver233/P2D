#include "full_cell_solver.h"
#include "../functions/functions.h"
#include <cstdlib>
#include <eigen3/Eigen/Sparse>
#include <cstdio>
#include <iostream>

inline double clamp(double x, double lower, double upper) {
    return x < lower ? lower : (x > upper ? upper : x);
}

/*
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
*/

void full_cell_solver::calc(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> c_s, double time, bool do_print) {
    int iter_time = 0;
    double first_norm, first_delta_norm;
    double res_norm = 1.0;
    double rel_tol = 1.0;
    double rel_delta = 1.0;
    MatrixXd du = MatrixXd::Zero(dof.dof_cnt, 1);
    std::vector<Eigen::Triplet<double> > coeff;
    //std::vector<double> temp;
    this->current_time = time;

    anode_particle.pre_calc(c_s);
    cathode_particle.pre_calc(c_s);
    //printf("Step\tIter\tRelTol\tDelta\n");

    while (iter_time < iter && rel_delta > tolerance) {
        Eigen::SparseMatrix<double> k(dof.dof_cnt, dof.dof_cnt);
        VectorXd res(dof.dof_cnt);
        coeff.clear();
        temp.clear();

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
        MatrixXd delta = -solver.solve(res);
        du += delta;
        u += delta;

        if (step != 0) {
            if (settings::calc_temperature || settings::use_adaptive_time_step) {
                double j1, j2;
                j1 = -anode_particle.calc<false>(c_s, u, 1, mesh, dof);
                j2 = -cathode_particle.calc<false>(c_s, u, 2, mesh, dof);
                anode.dc_ssdj = j1;
                cathode.dc_ssdj = j2;
            } else {
                anode_particle.calc<true>(c_s, u, 1, mesh, dof);
                cathode_particle.calc<true>(c_s, u, 2, mesh, dof);
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
        printf("%-8d%-8d%1.5lf %1.5lf %e\n", step, iter_time, rel_tol, rel_delta, res_norm);

        iter_time++;
    }
    if (do_print) this->print_detail();
    step++;
}


void full_cell_solver::print_detail() {
    /*
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
    Q_ohm.append(q_ohm, this->current_time);
    Q_rxn.append(q_rxn, this->current_time);
    Q_rev.append(q_rev, this->current_time);
    eta.append(v_eta, this->current_time);
    */
}
