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
    auto s = full_cell_solver(an, ca, ancoll, cacoll, coord, particle_coord);
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
            s.calc(u, c_s, i * constant::dt, false);
            if (i * constant::dt - s.next_detail_time > -0.001) {
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
                s.next_detail_time += constant::output_interval;
            }
        }
    } else {
        double current_time = 0.0;
        double max_dt = constant::output_interval;
        double cur_dt = max_dt;
        const int max_depth = 100;
        const double thres = 0.1;
        VectorXd mock_u;
        MatrixXd mock_c_s;
        VectorXd mock_u_2;
        MatrixXd mock_c_s_2;

        s.calc(u, c_s, 0, false);
        phi_l.append(u(Eigen::seq(0, pt_size - 1)), current_time);
        c_l.append(u(Eigen::seq(pt_size, 2 * pt_size - 1)), current_time);
        phi_s_negative.append(u(Eigen::seq(2 * pt_size, 2 * pt_size + (an - ancoll))), current_time);
        phi_s_positive.append(u(Eigen::seq(2 * pt_size + (an - ancoll) + 1, 2 * pt_size + eff_size - 1)),
                              current_time);
        j_negative.append(u(Eigen::seq(2 * pt_size + eff_size, 2 * pt_size + eff_size + (an - ancoll)))
                          * constant::F * 3 * constant::epsilon_s_an / constant::r_p * constant::j_ref,
                          current_time);
        j_positive.append(
            u(Eigen::seq(2 * pt_size + eff_size + (an - ancoll) + 1, 2 * pt_size + 2 * eff_size - 1))
            * constant::F * 3 * constant::epsilon_s_ca / constant::r_p * constant::j_ref, current_time);
        css_negative.append(c_s(
                                Eigen::seq(constant::particle_segment,
                                           (an - ancoll + 1) * (constant::particle_segment + 1),
                                           constant::particle_segment + 1), 0), current_time);
        css_positive.append(c_s(
                                Eigen::seq(
                                    (an - ancoll + 1) * (constant::particle_segment + 1) +
                                    constant::particle_segment,
                                    eff_size * (constant::particle_segment + 1),
                                    constant::particle_segment + 1), 0), current_time);
        if (settings::calc_temperature)
            temp.append(MatrixXd::Ones(1, 1) * u(2 * pt_size + 2 * eff_size), current_time);
        voltage.append(MatrixXd::Ones(1, 1) * (u(2 * pt_size + eff_size - 1) - u(2 * pt_size)),
                       current_time);
        s.print_detail();
        s.next_detail_time += constant::output_interval;
        while (current_time < constant::finish_time) {
            double diff_norm = 1e10;
            int depth = 0;
            double supposed_dt = std::min(cur_dt, std::min(max_dt, s.next_detail_time - current_time));
            double actual_dt = 0.0;
            while (diff_norm >= thres && depth <= max_depth) {
                mock_u = u;
                mock_c_s = c_s;
                mock_u_2 = u;
                mock_c_s_2 = c_s;
                actual_dt = supposed_dt;

                constant::dt = actual_dt;
                s.calc(mock_u, mock_c_s, current_time + actual_dt, false);

                constant::dt = actual_dt / 2;
                s.calc(mock_u_2, mock_c_s_2, current_time + actual_dt / 2, false);
                s.calc(mock_u_2, mock_c_s_2, current_time + actual_dt, false);

                diff_norm = (mock_u - mock_u_2).norm();
                std::cout << diff_norm << std::endl;
                depth++;
                supposed_dt /= 2.0;
            }
            if (depth > max_depth) {
                std::cerr << "Convergence failed" << std::endl;
                exit(1);
            } else if (depth != 1) {
                cur_dt = actual_dt;
            } else {
                cur_dt = std::min(max_dt, actual_dt * 2);
            }
            u = mock_u_2;
            c_s = mock_c_s_2;
            current_time += actual_dt;
            if (current_time - s.next_detail_time > -0.001) {
                phi_l.append(u(Eigen::seq(0, pt_size - 1)), current_time);
                c_l.append(u(Eigen::seq(pt_size, 2 * pt_size - 1)), current_time);
                phi_s_negative.append(u(Eigen::seq(2 * pt_size, 2 * pt_size + (an - ancoll))), current_time);
                phi_s_positive.append(u(Eigen::seq(2 * pt_size + (an - ancoll) + 1, 2 * pt_size + eff_size - 1)),
                                      current_time);
                j_negative.append(u(Eigen::seq(2 * pt_size + eff_size, 2 * pt_size + eff_size + (an - ancoll)))
                                  * constant::F * 3 * constant::epsilon_s_an / constant::r_p * constant::j_ref,
                                  current_time);
                j_positive.append(
                    u(Eigen::seq(2 * pt_size + eff_size + (an - ancoll) + 1, 2 * pt_size + 2 * eff_size - 1))
                    * constant::F * 3 * constant::epsilon_s_ca / constant::r_p * constant::j_ref, current_time);
                css_negative.append(c_s(
                                        Eigen::seq(constant::particle_segment,
                                                   (an - ancoll + 1) * (constant::particle_segment + 1),
                                                   constant::particle_segment + 1), 0), current_time);
                css_positive.append(c_s(
                                        Eigen::seq(
                                            (an - ancoll + 1) * (constant::particle_segment + 1) +
                                            constant::particle_segment,
                                            eff_size * (constant::particle_segment + 1),
                                            constant::particle_segment + 1), 0), current_time);
                if (settings::calc_temperature)
                    temp.append(MatrixXd::Ones(1, 1) * u(2 * pt_size + 2 * eff_size), current_time);
                voltage.append(MatrixXd::Ones(1, 1) * (u(2 * pt_size + eff_size - 1) - u(2 * pt_size)),
                               current_time);
                s.print_detail();
                s.next_detail_time += constant::output_interval;
            }
            std::cout << cur_dt << std::endl;
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
