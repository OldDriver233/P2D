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

    for(int i = 0; i < pt_size; i++) {
        u(i) = -uoc<1>(constant::c_int_an / constant::c_max_an);
    }
    //#pragma omp parallel for
    for(int i = pt_size; i < 2 * pt_size; i++) {
        u(i) = 1;
    }
    //#pragma omp parallel for
    for(int i = 2 * pt_size; i < 2 * pt_size + an - ancoll + 1; i++) {
        u(i) = 0;
    }
    //#pragma omp parallel for
    for(int i = 2 * pt_size + an - ancoll + 1; i < 2 * pt_size + eff_size; i++) {
        u(i) = uoc<2>(constant::c_int_ca / constant::c_max_ca) - uoc<1>(constant::c_int_an / constant::c_max_an);
    }
    if (settings::calc_temperature) {
        for(int i = 2 * pt_size + 2 * eff_size; i < 2 * pt_size + 2 * eff_size + all_size; i++) {
            u(i) = constant::t_ref;
        }
    }
    for(int i = 0; i < (an - ancoll + 1) * (constant::particle_segment + 1); i++) {
        c_s(i) = constant::c_int_an / constant::c_max_an;
    }
    for(int i = (an - ancoll + 1) * (constant::particle_segment + 1); i < eff_size * (constant::particle_segment + 1); i++) {
        c_s(i) = constant::c_int_ca / constant::c_max_ca;
    }

    std::vector<double> delta_u;
    std::vector<double> c_star;
    std::vector<double> voltage;
    std::vector<double> temperature;
    output_manager phi_l(coord(Eigen::seq(ancoll, cacoll)));
    output_manager c_l(coord(Eigen::seq(ancoll, cacoll)));
    output_manager phi_s_negative(coord(Eigen::seq(ancoll, an)));
    output_manager phi_s_positive(coord(Eigen::seq(ca, cacoll)));
    output_manager j_negative(coord(Eigen::seq(ancoll, an)));
    output_manager j_positive(coord(Eigen::seq(ca, cacoll)));
    output_manager css_negative(coord(Eigen::seq(ancoll, an)));
    output_manager css_positive(coord(Eigen::seq(ca, cacoll)));
    output_manager temp(MatrixXd::Zero(1, 1));
    for(int i = 0; i <= constant::step; i++) {
        s.calc(u, c_s);
        //std::cout<<"Step: "<<i<<std::endl;
        delta_u.push_back(u(pt_size - 1) - u(0));
        voltage.push_back(u(2 * pt_size + eff_size - 1) - u(2 * pt_size));
        if (settings::calc_temperature)
            temperature.push_back((u(2 * pt_size + 2 * eff_size) + u(2 * pt_size + 2 * all_size - 1)) / 2);
        if (static_cast<int>(i * constant::dt) % 36 == 0) {
            phi_l.append(u(Eigen::seq(0, pt_size - 1)), constant::dt * i);
            c_l.append(u(Eigen::seq(pt_size, 2 * pt_size - 1)), constant::dt * i);
            phi_s_negative.append(u(Eigen::seq(2 * pt_size, 2 * pt_size + (an - ancoll))), constant::dt * i);
            phi_s_positive.append(u(Eigen::seq(2 * pt_size + (an - ancoll) + 1, 2 * pt_size + eff_size - 1)), constant::dt * i);
            j_negative.append(u(Eigen::seq(2 * pt_size + eff_size, 2 * pt_size + eff_size + (an - ancoll)))
                * constant::F * 3 * constant::epsilon_s_an / constant::r_p * constant::j_ref, constant::dt * i);
            j_positive.append(u(Eigen::seq(2 * pt_size + eff_size + (an - ancoll) + 1, 2 * pt_size + 2 * eff_size - 1))
                * constant::F * 3 * constant::epsilon_s_ca / constant::r_p * constant::j_ref, constant::dt * i);
            css_negative.append(c_s(
                Eigen::seq(constant::particle_segment,
                    (an - ancoll + 1) * (constant::particle_segment + 1),
                    constant::particle_segment + 1), 0), constant::dt * i);
            css_positive.append(c_s(
                Eigen::seq((an - ancoll + 1) * (constant::particle_segment + 1) + constant::particle_segment,
                    eff_size * (constant::particle_segment + 1),
                    constant::particle_segment + 1), 0), constant::dt * i);
            if (settings::calc_temperature)
                temp.append(MatrixXd::Ones(1, 1) * u(2 * pt_size + 2 * eff_size), constant::dt * i);
        }

    }
    phi_l.write_to_csv("output/phi_l.csv");
    c_l.write_to_csv("output/c_l.csv");
    phi_s_negative.write_to_csv("output/phi_neg.csv");
    phi_s_positive.write_to_csv("output/phi_pos.csv");
    temp.write_to_csv("output/temp.csv");
    j_positive.write_to_csv("output/j_pos.csv");
    j_negative.write_to_csv("output/j_neg.csv");
    s.Q_ohm.write_to_csv("output/q_ohm.csv");
    s.Q_rxn.write_to_csv("output/q_rxn.csv");
    s.Q_rev.write_to_csv("output/q_rev.csv");
    s.eta.write_to_csv("output/eta.csv");
    css_negative.write_to_csv("output/css_neg.csv");
    css_positive.write_to_csv("output/css_pos.csv");

    /*
    std::cout<<u<<std::endl;
    if (settings::calc_temperature) {
        for(int i = 0; i < all_size; i++) {
            std::cout<<u(2 * eff_size + 2 * pt_size + i) - constant::t_ref<<std::endl;
        }
    }
    std::cout<<std::endl;
    for(int i = 0; i < eff_size; i++) {
        std::cout<<c_s((i + 1) * (constant::particle_segment + 1) - 1)<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<c_s.block(0, 0, constant::particle_segment + 1, 1)<<std::endl;
    */

    /*
    std::cout<<"Writing to redis"<<std::endl;
    auto redis = redis_connector();
    redis.del("voltage");
    redis.del("delta_u");
    redis.del("u");
    redis.del("c_star");
    redis.del("temp");
    redis.del("phi_e");
    for(int i = 0; i < voltage.size(); i++) {
        redis.rpush("voltage", std::to_string(voltage[i]));
    }
    for(int i = 0; i < delta_u.size(); i++) {
        redis.rpush("delta_u", std::to_string(delta_u[i]));
    }
    for(int i = 0; i < c_star.size(); i++) {
        redis.rpush("c_star", std::to_string(c_star[i] * constant::c_max_ca / 10000));
    }
    for(int i = 0; i < pt_size; i++) {
        redis.rpush("u", std::to_string(u(i + pt_size) * constant::ce_int));
    }
    for (int i = 0; i < pt_size; i++) {
        redis.rpush("phi_e", std::to_string(u(i)));
    }
    if (settings::calc_temperature) {
        for(int i = 0; i < temperature.size(); i++) {
            redis.rpush("temp", std::to_string(temperature[i]));
        }
    }
    redis.set("last_update_at", std::to_string(get_timestamp()));
    */
}

/*
void test_particle() {
    constant::read();
    VectorXd coord = VectorXd::LinSpaced(constant::particle_segment + 1, 0.0, 1.0);
    auto p = particle_solver(coord, constant::ds_an, constant::c_max_an);
    std::cout<<p.assembled_A<<std::endl;
    std::cout<<p.point_coord<<std::endl;
    //std::cout<<p.inv_AB<<std::endl;
    std::cout<<p.j_coeff<<std::endl;
}
*/

int main() {
    calc_cell();
    //test_particle();
}
