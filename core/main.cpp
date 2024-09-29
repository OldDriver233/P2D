#include <fstream>
#include <iostream>
#include "constants/constant.h"
#include "functions/functions.h"
#include "mesh/mesh.h"
#include "mesh/mesh_reader.h"
#include "particle/particle_solver.h"
#include "full_cell/full_cell_solver.h"
#include "io/coord_reader.h"
#include "io/redis_connector.h"
#include <eigen3/Eigen/Dense>
#include <vector>
#include <chrono>

using Eigen::MatrixXd;
using Eigen::VectorXd;

//using namespace sw::redis;

void check_mesh() {
    mesh m;
    auto mr = mesh_reader("macro2d.msh");
    mr.read(m);
}

std::time_t get_timestamp() {
    auto tp = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
    auto stamp = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch());
    std::time_t timestamp = stamp.count();
    return timestamp;
}

void calc_cell() {
    auto [coord, an, ca] = coord_reader("coord.json");
    int pt_size = coord.size();
    int eff_size = pt_size - (ca - an - 1);
    constant::read();
    auto s = full_cell_solver(an, ca, coord);
    MatrixXd u = MatrixXd::Zero(2 * pt_size + 4 * eff_size, 1);
    for(int i = pt_size; i < 2 * pt_size; i++) {
        u(i) = 1;
    }
    for(int i = 2 * pt_size; i < 2 * pt_size + an + 1; i++) {
        u(i) = uoc(constant::c_int_an / constant::c_max_an, 1);
    }
    for(int i = 2 * pt_size + an + 1; i < 2 * pt_size + eff_size; i++) {
        u(i) = uoc(constant::c_int_ca / constant::c_max_ca, 2);
    }
    for(int i = 2 * pt_size + 2 * eff_size; i < 2 * pt_size + 2 * eff_size + an + 1; i++) {
        u(i) = constant::c_int_an / constant::c_max_an;
    }
    for(int i = 2 * pt_size + 2 * eff_size + an + 1; i < 2 * pt_size + 3 * eff_size; i++) {
        u(i) = constant::c_int_ca / constant::c_max_ca;
    }
    for(int i = 2 * pt_size + 3 * eff_size; i < 2 * pt_size + 3 * eff_size + an + 1; i++) {
        u(i) = constant::c_int_an / constant::c_max_an;
    }
    for(int i = 2 * pt_size + 3 * eff_size + an + 1; i < 2 * pt_size + 4 * eff_size; i++) {
        u(i) = constant::c_int_ca / constant::c_max_ca;
    }

    std::vector<double> delta_u;
    std::vector<double> c_star;
    std::vector<double> voltage;
    for(int i = 0; i <= constant::step; i++) {
        s.calc(u);
        //std::cout<<"Step: "<<i<<std::endl;
        delta_u.push_back(u(pt_size - 1) - u(0));
        c_star.push_back(u(2 * pt_size + 4 * eff_size - 1));
        voltage.push_back(u(2 * pt_size + eff_size - 1) - u(2 * pt_size));
    }
    std::cout<<u<<std::endl;

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


void test() {
    int pt_size = 16;
    int an = -1, ca = -1;
    constant::read();
    VectorXd coord = VectorXd::LinSpaced(pt_size, 0, 1);
    auto s = full_cell_solver(an, ca, VectorXd::LinSpaced(pt_size, 0, 1));
    MatrixXd u = MatrixXd::Zero(2 * pt_size, 1);
    for(int i = pt_size; i < 2 * pt_size; i++) {
        u(i) = 1000.0;
    }
    for(int i = 0; i < constant::step; i++) {
        s.calc(u);
    }
    std::cout<<u<<std::endl;
}

int main() {
    //calc_separator();
    calc_cell();
    //test();
}
