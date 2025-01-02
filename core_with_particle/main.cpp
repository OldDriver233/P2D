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
    MatrixXd u = MatrixXd::Zero(2 * pt_size + 2 * eff_size + all_size, 1);
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
    for(int i = 2 * pt_size + 2 * eff_size; i < 2 * pt_size + 2 * eff_size + all_size; i++) {
        u(i) = 297.0;
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
    for(int i = 0; i <= constant::step; i++) {
        s.calc(u, c_s);
        //std::cout<<"Step: "<<i<<std::endl;
        delta_u.push_back(u(pt_size - 1) - u(0));
        voltage.push_back(u(2 * pt_size + eff_size - 1) - u(2 * pt_size));
    }
    std::cout<<u<<std::endl;
    for(int i = 0; i < all_size; i++) {
        std::cout<<u(2 * eff_size + 2 * pt_size + i) - 297<<std::endl;
    }
    std::cout<<std::endl;
    for(int i = 0; i < eff_size; i++) {
        std::cout<<c_s((i + 1) * (constant::particle_segment + 1) - 1)<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<c_s.block(0, 0, constant::particle_segment + 1, 1)<<std::endl;

    std::cout<<"Writing to redis"<<std::endl;
    auto redis = redis_connector();
    redis.del("voltage");
    redis.del("delta_u");
    redis.del("u");
    redis.del("c_star");
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
    redis.set("last_update_at", std::to_string(get_timestamp()));
}

void test_particle() {
    constant::read();
    VectorXd coord = VectorXd::LinSpaced(constant::particle_segment + 1, 0.0, 1.0);
    auto p = particle_solver(coord, constant::ds_an, constant::c_max_an);
    std::cout<<p.assembled_A<<std::endl;
    std::cout<<p.point_coord<<std::endl;
    //std::cout<<p.inv_AB<<std::endl;
    std::cout<<p.j_coeff<<std::endl;
}

int main() {
    calc_cell();
    //test_particle();
}
