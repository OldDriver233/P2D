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

void full_cell_solver::calc(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> c_s, std::vector<double>& v_stress, double time, bool do_print) {
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
    Eigen::SparseMatrix<double> k(dof.dof_cnt, dof.dof_cnt);
    if (step == 0 && settings::stress_analysis) {
        stress.generate(stress_mat_coeff, local_stress_mat_coeff);
        stress_mat.setFromTriplets(local_stress_mat_coeff.begin(), local_stress_mat_coeff.end());
    }

    while (iter_time < iter && rel_delta > tolerance) {
        k.setZero();
        VectorXd res = VectorXd::Zero(dof.dof_cnt);
        coeff.clear();
        temp.clear();
        std::vector<double> avg_conc(dof.particle_to_node.size());
        anode_particle.get_average_concentration(c_s, avg_conc, 1, mesh, dof);
        cathode_particle.get_average_concentration(c_s, avg_conc, 2, mesh, dof);

        if (settings::calc_temperature) {
            this->anode.generate<true>(u, du, c_s, v_stress, avg_conc, coeff, res, step == 0, temp);
            this->sep.generate<true>(u, du, c_s, coeff, res, step == 0, temp);
            this->cathode.generate<true>(u, du, c_s, v_stress, avg_conc, coeff, res, step == 0, temp);
            this->anode_collector.generate(u, du, c_s, coeff, res, step == 0);
            this->cathode_collector.generate(u, du, c_s, coeff, res, step == 0);
        } else {
            this->anode.generate<false>(u, du, c_s, v_stress, avg_conc, coeff, res, step == 0, temp);
            this->sep.generate<false>(u, du, c_s, coeff, res, step == 0, temp);
            this->cathode.generate<false>(u, du, c_s, v_stress, avg_conc, coeff, res, step == 0, temp);
        }
        if (settings::stress_analysis) coeff.insert(coeff.end(), stress_mat_coeff.begin(), stress_mat_coeff.end());
        k.setFromTriplets(coeff.begin(), coeff.end());
        if (settings::stress_analysis) stress.generate_residue(u, avg_conc, res, stress_mat);
        apply_boundary(u, k, res, false);

        VectorXd delta;
        if (step != 0 || true) {
            solver.compute(k);
            delta = -solver.solve(res);
            std::cout<<std::setw(12);
        } else {
            /*
            solver.compute(k);
            VectorXd direction = -solver.solve(res);
            double alpha = 0.01, prev_alpha = 1;
            int inner_iter = 0;

            while (abs(alpha - prev_alpha) > 0.005 && inner_iter < 50) {
                prev_alpha = alpha;
                k.setZero();
                res = VectorXd::Zero(dof.dof_cnt);
                coeff.clear();
                VectorXd u_a = u + alpha * direction;
                VectorXd du_a = du + alpha * direction;

                if (settings::calc_temperature) {
                    this->anode.generate<true>(u_a, du_a, c_s, coeff, res, step == 0, temp);
                    this->sep.generate<true>(u_a, du_a, c_s, coeff, res, step == 0, temp);
                    this->cathode.generate<true>(u_a, du_a, c_s, coeff, res, step == 0, temp);
                    this->anode_collector.generate(u_a, du_a, c_s, coeff, res, step == 0);
                    this->cathode_collector.generate(u_a, du_a, c_s, coeff, res, step == 0);
                } else {
                    this->anode.generate<false>(u_a, du_a, c_s, coeff, res, step == 0, temp);
                    this->sep.generate<false>(u_a, du_a, c_s, coeff, res, step == 0, temp);
                    this->cathode.generate<false>(u_a, du_a, c_s, coeff, res, step == 0, temp);
                }
                k.setFromTriplets(coeff.begin(), coeff.end());
                apply_boundary(u_a, k, res, false);

                alpha = alpha - direction.dot(res) / (direction.transpose().dot(k * direction));

                std::cout<<alpha<<" "<<prev_alpha<<" "<<res.norm()<<std::endl;
                inner_iter++;
            }

            alpha = std::clamp(alpha, .001, 1.0);
            delta = direction * alpha;
            */
        }
        du += delta;
        u += delta;
        //std::cout<<res.transpose()<<std::endl;
        //std::cout<<delta.transpose()<<std::endl;
        //std::cout<<u.transpose()<<std::endl;
        //std::cout<<std::endl;

        if (step != 0) {
            if (settings::calc_temperature || settings::use_adaptive_time_step || settings::stress_analysis) {
                double j1, j2;
                if (!settings::stress_analysis) {
                    j1 = -anode_particle.calc(c_s, u, 1, mesh, dof);
                    j2 = -cathode_particle.calc(c_s, u, 2, mesh, dof);
                } else {
                    j1 = -anode_particle.calc_stress(c_s, u, 1, mesh, dof);
                    j2 = -cathode_particle.calc_stress(c_s, u, 2, mesh, dof);
                }
                anode.dc_ssdj = j1;
                cathode.dc_ssdj = j2;
            } else {
                anode_particle.calc(c_s, u, 1, mesh, dof);
                cathode_particle.calc(c_s, u, 2, mesh, dof);
            }
        }
        if (settings::stress_analysis) stress.stress_output(u, v_stress);
        double norm = delta.norm();
        res_norm = res.norm();
        if (iter_time <= 1) {
            first_norm = res_norm;
            first_delta_norm = norm;
        } else {
            rel_tol = res_norm / first_norm;
            rel_delta = norm / first_delta_norm;
        }
        if (iter_time == 0) printf("%-8d%-8d-       -       %e %e\n", step, iter_time, res_norm, norm);
        else printf("%-8d%-8d%1.5lf %1.5lf %e %e\n", step, iter_time, rel_tol, rel_delta, res_norm, norm);

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
