#include <fstream>
#include <iostream>
#include "constants/constant.h"
#include "mesh/mesh.h"
#include "mesh/mesh_reader.h"
#include "particle/particle_solver.h"
#include <eigen3/Eigen/Dense>
#include <sw/redis++/redis.h>
#include <vector>
#include <sw/redis++/redis++.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace sw::redis;

void check_mesh() {
    mesh m;
    auto mr = mesh_reader("macro2d.msh");
    mr.read(m);
}

int main()
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
    pipe.exec();
    return 0;
}

/*
namespace py = pybind11;

PYBIND11_MODULE(fem, m) {
    py::class_<particle_solver>(m, "ParticleSolver")
            .def(py::init<const VectorXd&>())
            .def("calc", &particle_solver::calc);

}
 */