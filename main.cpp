#include <fstream>
#include <iostream>
#include "mesh/mesh.h"
#include "mesh/mesh_reader.h"
#include "particle/particle_solver.h"
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

void check_mesh() {
    mesh m;
    auto mr = mesh_reader("macro2d.msh");
    mr.read(m);
}

int main()
{
    if(true) {
        check_mesh();
        return 0;
    }
    VectorXd coord = VectorXd::LinSpaced(101, 0, 1);
    auto s = particle_solver(VectorXd::LinSpaced(101, 0, 1));
    MatrixXd u = MatrixXd::Zero(101, 1);
    for(int i = 0; i < 10; i++) {
        s.calc(u);
    }

    std::fstream fs;
    fs.open("out.csv", std::ios::out);
    fs<<"x, y"<<std::endl;
    for(int i = 0; i < 101; i++) {
        fs<<coord(i)<<", "<<u(i)<<std::endl;
    }
    fs.close();
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