#include <fstream>
#include <iostream>
#include "constants/constant.h"
#include "functions/functions.h"
#include "mesh/mesh.h"
#include "mesh/mesh_reader.h"
#include "particle/particle_solver.h"
#include "full_cell/full_cell_solver.h"
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

/*
void calc_particle()
{
    //if(true) {
    //    check_mesh();
    //    return 0;
    //}
    std::vector<double> us;
    constant::read();
    VectorXd coord = VectorXd::LinSpaced(101, 0, 1);
    auto s = particle_solver(VectorXd::LinSpaced(101, 0, 1));
    MatrixXd u = MatrixXd::Ones(101, 1) * constant::initial_c / constant::max_c;
    for(int i = 0; i < constant::step; i++) {
        us.push_back(u(100));
        s.calc(u);
    }
    us.push_back(u(100));

    std::cout<<"Calculation complete. Writing to redis..."<<std::endl;
    
    auto redis = Redis("tcp://127.0.0.1:6379");
    auto pipe = redis.pipeline();
    pipe.del("us");
    for(int i = 0; i < us.size(); i++) {
        pipe.rpush("us", std::to_string(us[i]));
    }
    pipe.set("us:last_update_at", std::to_string(get_timestamp()));
    pipe.exec();
}
*/
void calc_separator() {
    
}


void calc_cell() {
    int pt_size = 46;
    int an = 20, ca = 25;
    int eff_size = pt_size - (ca - an - 1);
    constant::read();
    VectorXd coord = VectorXd::LinSpaced(pt_size, 0, 225);
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
    for(int i = 0; i < constant::step; i++) {
        s.calc(u);
        std::cout<<"Step: "<<i<<std::endl;
        delta_u.push_back(u(pt_size - 1) - u(0));
    }
    std::cout<<u<<std::endl;

    std::cout<<"Writing to redis"<<std::endl;
    auto redis = redis_connector();
    redis.del("delta_u");
    redis.del("u");
    for(int i = 0; i < delta_u.size(); i++) {
        redis.rpush("delta_u", std::to_string(delta_u[i]));
    }
    for(int i = 0; i < pt_size; i++) {
        redis.rpush("u", std::to_string(u(i + pt_size) - 1));
    }
    redis.set("last_update_at", std::to_string(get_timestamp()));
    /*
    auto redis = Redis("tcp://127.0.0.1:6379");
    auto pipe = redis.pipeline();
    pipe.del("us");
    for(int i = 0; i < u.size() / 2; i++) {
        pipe.rpush("us", std::to_string(u(i)));
    }
    pipe.set("us:last_update_at", std::to_string(get_timestamp()));
    pipe.exec();
    */
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

/*
namespace py = pybind11;

PYBIND11_MODULE(fem, m) {
    py::class_<particle_solver>(m, "ParticleSolver")
            .def(py::init<const VectorXd&>())
            .def("calc", &particle_solver::calc);

}
 */