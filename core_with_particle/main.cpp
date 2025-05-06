#include <iostream>
#include "constants/constant.h"
#include "functions/functions.h"
#include "full_cell/full_cell_solver.h"
#include "full_cell/particle/particle_solver.h"
#include "io/coord_reader.h"
#include "io/redis_connector.h"
#include "io/settings/settings.h"
#include <eigen3/Eigen/Dense>
#include <ostream>
#include <vector>
#include <chrono>

#include "io/output/output_manager.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

//using namespace sw::redis;

std::time_t get_timestamp() {
    auto tp = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
    auto stamp = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch());
    std::time_t timestamp = stamp.count();
    return timestamp;
}

void calc_cell() {
    settings::read("settings.json");
    auto [coord, an, ca, ancoll, cacoll] = coord_reader();
    int pt_size = cacoll - ancoll + 1;
    int eff_size = pt_size - (ca - an - 1);
    int all_size = coord.size();
    constant::read();
    MatrixXd particle_coord = VectorXd::LinSpaced(constant::particle_segment + 1, 0.0, 1.0);
    StepControl step_control;
    auto s = full_cell_solver(an, ca, ancoll, cacoll, coord, particle_coord, &step_control);
    VectorXd u;
    if (settings::calc_temperature) u = MatrixXd::Zero(2 * pt_size + 2 * eff_size + all_size, 1);
    else u = MatrixXd::Zero(2 * pt_size + 2 * eff_size, 1);
    MatrixXd c_s = MatrixXd::Zero(eff_size * (constant::particle_segment + 1), 1);

    for (int i = 0; i < pt_size; i++) {
        u(i) = -uoc<1>(constant::c_int_an / constant::c_max_an);
    }
    //#pragma omp parallel for
    for (int i = pt_size; i < 2 * pt_size; i++) {
        u(i) = 1;
    }
    //#pragma omp parallel for
    for (int i = 2 * pt_size; i < 2 * pt_size + an - ancoll + 1; i++) {
        u(i) = 0;
    }
    //#pragma omp parallel for
    for (int i = 2 * pt_size + an - ancoll + 1; i < 2 * pt_size + eff_size; i++) {
        u(i) = uoc<2>(constant::c_int_ca / constant::c_max_ca) - uoc<1>(constant::c_int_an / constant::c_max_an);
    }
    if (settings::calc_temperature) {
        for (int i = 2 * pt_size + 2 * eff_size; i < 2 * pt_size + 2 * eff_size + all_size; i++) {
            u(i) = constant::t_ref;
        }
    }
    for (int i = 0; i < (an - ancoll + 1) * (constant::particle_segment + 1); i++) {
        c_s(i) = constant::c_int_an / constant::c_max_an;
    }
    for (int i = (an - ancoll + 1) * (constant::particle_segment + 1); i < eff_size * (constant::particle_segment + 1);
         i++) {
        c_s(i) = constant::c_int_ca / constant::c_max_ca;
         }

    output_manager phi_l(coord(Eigen::seq(ancoll, cacoll)));
    output_manager c_l(coord(Eigen::seq(ancoll, cacoll)));
    output_manager phi_s_negative(coord(Eigen::seq(ancoll, an)));
    output_manager phi_s_positive(coord(Eigen::seq(ca, cacoll)));
    output_manager j_negative(coord(Eigen::seq(ancoll, an)));
    output_manager j_positive(coord(Eigen::seq(ca, cacoll)));
    output_manager css_negative(coord(Eigen::seq(ancoll, an)));
    output_manager css_positive(coord(Eigen::seq(ca, cacoll)));
    output_manager temp(MatrixXd::Zero(1, 1));
    output_manager voltage(MatrixXd::Zero(1, 1));
    if (!settings::use_adaptive_time_step) {
        for (int i = 0; i <= constant::step; i++) {
            step_control.dt_now = constant::dt;
            s.calc(u, c_s, i * constant::dt, false);
            if (i * constant::dt - s.step_control->next_output > -0.001) {
                phi_l.append(u(Eigen::seq(0, pt_size - 1)), constant::dt * i);
                c_l.append(u(Eigen::seq(pt_size, 2 * pt_size - 1)), constant::dt * i);
                phi_s_negative.append(u(Eigen::seq(2 * pt_size, 2 * pt_size + (an - ancoll))), constant::dt * i);
                phi_s_positive.append(u(Eigen::seq(2 * pt_size + (an - ancoll) + 1, 2 * pt_size + eff_size - 1)),
                                      constant::dt * i);
                j_negative.append(u(Eigen::seq(2 * pt_size + eff_size, 2 * pt_size + eff_size + (an - ancoll)))
                                  * constant::F * 3 * constant::epsilon_s_an / constant::r_p * constant::j_ref,
                                  constant::dt * i);
                j_positive.append(
                    u(Eigen::seq(2 * pt_size + eff_size + (an - ancoll) + 1, 2 * pt_size + 2 * eff_size - 1))
                    * constant::F * 3 * constant::epsilon_s_ca / constant::r_p * constant::j_ref, constant::dt * i);
                css_negative.append(c_s(
                                        Eigen::seq(constant::particle_segment,
                                                   (an - ancoll + 1) * (constant::particle_segment + 1),
                                                   constant::particle_segment + 1), 0), constant::dt * i);
                css_positive.append(c_s(
                                        Eigen::seq(
                                            (an - ancoll + 1) * (constant::particle_segment + 1) +
                                            constant::particle_segment,
                                            eff_size * (constant::particle_segment + 1),
                                            constant::particle_segment + 1), 0), constant::dt * i);
                if (settings::calc_temperature)
                    temp.append(MatrixXd::Ones(1, 1) * u(2 * pt_size + 2 * eff_size), constant::dt * i);
                voltage.append(MatrixXd::Ones(1, 1) * (u(2 * pt_size + eff_size - 1) - u(2 * pt_size)),
                               constant::dt * i);
                s.print_detail();
                s.step_control->next_output += constant::output_interval;
            }
        }
    } else {
        s.calc(u, c_s, 0, false);
        phi_l.append(u(Eigen::seq(0, pt_size - 1)), 0);
        c_l.append(u(Eigen::seq(pt_size, 2 * pt_size - 1)), 0);
        phi_s_negative.append(u(Eigen::seq(2 * pt_size, 2 * pt_size + (an - ancoll))), 0);
        phi_s_positive.append(u(Eigen::seq(2 * pt_size + (an - ancoll) + 1, 2 * pt_size + eff_size - 1)),
                              0);
        j_negative.append(u(Eigen::seq(2 * pt_size + eff_size, 2 * pt_size + eff_size + (an - ancoll)))
                          * constant::F * 3 * constant::epsilon_s_an / constant::r_p * constant::j_ref,
                          0);
        j_positive.append(
            u(Eigen::seq(2 * pt_size + eff_size + (an - ancoll) + 1, 2 * pt_size + 2 * eff_size - 1))
            * constant::F * 3 * constant::epsilon_s_ca / constant::r_p * constant::j_ref, 0);
        css_negative.append(c_s(
                                Eigen::seq(constant::particle_segment,
                                           (an - ancoll + 1) * (constant::particle_segment + 1),
                                           constant::particle_segment + 1), 0), 0);
        css_positive.append(c_s(
                                Eigen::seq(
                                    (an - ancoll + 1) * (constant::particle_segment + 1) +
                                    constant::particle_segment,
                                    eff_size * (constant::particle_segment + 1),
                                    constant::particle_segment + 1), 0), 0);
        if (settings::calc_temperature)
            temp.append(MatrixXd::Ones(1, 1) * u(2 * pt_size + 2 * eff_size), 0);
        voltage.append(MatrixXd::Ones(1, 1) * (u(2 * pt_size + eff_size - 1) - u(2 * pt_size)),
                       0);
        s.print_detail();

        step_control.update_solution(u, c_s, 1);
        step_control.next_output = constant::output_interval;
        while (step_control.prev_time < constant::finish_time) {
            s.calc(u, c_s, step_control.prev_time + step_control.dt_now, false);
            step_control.update_solution(u, c_s, 0);
            step_control.update_dt();
            while (step_control.status == StepStatus::NeedRecalc) {
                u = step_control.u_hist[1];
                c_s = step_control.c_s_hist[1];
                s.calc(u, c_s, step_control.prev_time + step_control.dt_now, false);
                step_control.update_solution(u, c_s, 0);
                step_control.update_dt();
            }
            std::cout<<step_control.prev_time<<" "<<step_control.dt_proposed<<std::endl;
            if (step_control.status == StepStatus::ToNextOutput) {
                phi_l.append(u(Eigen::seq(0, pt_size - 1)), step_control.next_output);
                c_l.append(u(Eigen::seq(pt_size, 2 * pt_size - 1)), step_control.next_output);
                phi_s_negative.append(u(Eigen::seq(2 * pt_size, 2 * pt_size + (an - ancoll))), step_control.next_output);
                phi_s_positive.append(u(Eigen::seq(2 * pt_size + (an - ancoll) + 1, 2 * pt_size + eff_size - 1)),
                                      step_control.next_output);
                j_negative.append(u(Eigen::seq(2 * pt_size + eff_size, 2 * pt_size + eff_size + (an - ancoll)))
                                  * constant::F * 3 * constant::epsilon_s_an / constant::r_p * constant::j_ref,
                                  step_control.next_output);
                j_positive.append(
                    u(Eigen::seq(2 * pt_size + eff_size + (an - ancoll) + 1, 2 * pt_size + 2 * eff_size - 1))
                    * constant::F * 3 * constant::epsilon_s_ca / constant::r_p * constant::j_ref, step_control.next_output);
                css_negative.append(c_s(
                                        Eigen::seq(constant::particle_segment,
                                                   (an - ancoll + 1) * (constant::particle_segment + 1),
                                                   constant::particle_segment + 1), 0), step_control.next_output);
                css_positive.append(c_s(
                                        Eigen::seq(
                                            (an - ancoll + 1) * (constant::particle_segment + 1) +
                                            constant::particle_segment,
                                            eff_size * (constant::particle_segment + 1),
                                            constant::particle_segment + 1), 0), step_control.next_output);
                if (settings::calc_temperature)
                    temp.append(MatrixXd::Ones(1, 1) * u(2 * pt_size + 2 * eff_size), step_control.next_output);
                voltage.append(MatrixXd::Ones(1, 1) * (u(2 * pt_size + eff_size - 1) - u(2 * pt_size)),
                               step_control.next_output);
                s.print_detail();
                step_control.next_output += constant::output_interval;
            }
            step_control.step_forward();
        }
        if (step_control.status == StepStatus::ToNextOutput) {
                phi_l.append(u(Eigen::seq(0, pt_size - 1)), step_control.next_output);
                c_l.append(u(Eigen::seq(pt_size, 2 * pt_size - 1)), step_control.next_output);
                phi_s_negative.append(u(Eigen::seq(2 * pt_size, 2 * pt_size + (an - ancoll))), step_control.next_output);
                phi_s_positive.append(u(Eigen::seq(2 * pt_size + (an - ancoll) + 1, 2 * pt_size + eff_size - 1)),
                                      step_control.next_output);
                j_negative.append(u(Eigen::seq(2 * pt_size + eff_size, 2 * pt_size + eff_size + (an - ancoll)))
                                  * constant::F * 3 * constant::epsilon_s_an / constant::r_p * constant::j_ref,
                                  step_control.next_output);
                j_positive.append(
                    u(Eigen::seq(2 * pt_size + eff_size + (an - ancoll) + 1, 2 * pt_size + 2 * eff_size - 1))
                    * constant::F * 3 * constant::epsilon_s_ca / constant::r_p * constant::j_ref, step_control.next_output);
                css_negative.append(c_s(
                                        Eigen::seq(constant::particle_segment,
                                                   (an - ancoll + 1) * (constant::particle_segment + 1),
                                                   constant::particle_segment + 1), 0), step_control.next_output);
                css_positive.append(c_s(
                                        Eigen::seq(
                                            (an - ancoll + 1) * (constant::particle_segment + 1) +
                                            constant::particle_segment,
                                            eff_size * (constant::particle_segment + 1),
                                            constant::particle_segment + 1), 0), step_control.next_output);
                if (settings::calc_temperature)
                    temp.append(MatrixXd::Ones(1, 1) * u(2 * pt_size + 2 * eff_size), step_control.next_output);
                voltage.append(MatrixXd::Ones(1, 1) * (u(2 * pt_size + eff_size - 1) - u(2 * pt_size)),
                               step_control.next_output);
                s.print_detail();
            }
    }
    phi_l.write_to_csv("output/phi_l.csv");
    c_l.write_to_csv("output/c_l.csv");
    phi_s_negative.write_to_csv("output/phi_neg.csv");
    phi_s_positive.write_to_csv("output/phi_pos.csv");
    temp.write_to_csv("output/temp.csv");
    voltage.write_to_csv("output/voltage.csv");
    j_positive.write_to_csv("output/j_pos.csv");
    j_negative.write_to_csv("output/j_neg.csv");
    s.Q_ohm.write_to_csv("output/q_ohm.csv");
    s.Q_rxn.write_to_csv("output/q_rxn.csv");
    s.Q_rev.write_to_csv("output/q_rev.csv");
    s.eta.write_to_csv("output/eta.csv");
    css_negative.write_to_csv("output/css_neg.csv");
    css_positive.write_to_csv("output/css_pos.csv");
}

int main() {
    calc_cell();
    //test_particle();
}
