#include "particle_solver.h"
#include "stiffness_generator.h"
#include <cmath>
#include <eigen3/Eigen/src/Core/ArithmeticSequence.h>
#include <iostream>

double clamp(double x, double lower, double upper) {
    return x < lower ? lower : (x > upper ? upper : x);
}

double j_border(double c) {
    return (c - c * c) * 1e-7;
}

double d_j_border(double c) {
    return (1.0 - 2.0 * c) * 1e-7;
}

void particle_solver::calc(Eigen::Ref<MatrixXd> u) {
    long element_cnt = this->point_coord.size() - 1;
    result = VectorXd::Zero(element_cnt + 1);

    int iter_time = 0;
    double ratio = 999999.0;
    double first_norm;
    double eff = M_PI * 4.0 * r / d_ref;
    MatrixXd du = MatrixXd::Zero(element_cnt + 1, 1);

    while(iter_time < iter && ratio > tolerance) {
        this->gen.generate(u, du);
            
        gen.generated_K(element_cnt, element_cnt) -= d_j_border(u(element_cnt, 0)) * eff;
        gen.generated_res(element_cnt) -= j_border(u(element_cnt, 0)) * eff;

        MatrixXd delta = -gen.generated_K.fullPivLu().solve(gen.generated_res);
        du += delta;
        u += delta;
        double norm = delta.norm();
        if(iter_time != 0) {
            ratio = norm / first_norm;
        } else {
            first_norm = norm;
        }

        iter_time++;
    }
}
