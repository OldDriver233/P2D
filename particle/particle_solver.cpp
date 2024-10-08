#include "particle_solver.h"
#include "stiffness_generator.h"
#include <cmath>
#include <eigen3/Eigen/src/Core/ArithmeticSequence.h>
#include <iostream>


double clamp(double x, double lower, double upper) {
    return x < lower ? lower : (x > upper ? upper : x);
}

double ueq(double c, int type) {
    if(type == 1) return 0.194 + 1.5 * std::exp(-120.0 * c)
        + 0.0351 * std::tanh((c - 0.286) / 0.083)
        - 0.0045 * std::tanh((c - 0.849) / 0.119)
        - 0.035 * std::tanh((c - 0.9233) / 0.05)
        - 0.0147 * std::tanh((c - 0.5) / 0.034)
        - 0.102 * std::tanh((c - 0.194) / 0.142)
        - 0.022 * std::tanh((c - 0.9) / 0.0164)
        - 0.011 * std::tanh((c - 0.124) / 0.0226)
        + 0.0155 * std::tanh((c - 0.105) / 0.029);
    else return 2.16216 + 0.07645 * std::tanh(30.834 - 54.4806 * c)
        + 2.1581 * std::tanh(52.294 - 50.294 * c)
        - 0.14169 * std::tanh(11.0923 - 19.8543 * c)
        + 0.2051 * std::tanh(1.4684 - 5.4888 * c)
        + 0.2531 * std::tanh((-c + 0.56468) / 0.1316)
        - 0.02167 * std::tanh((c - 0.525) / 0.006);
}

double d_ueq(double c, int type) {
    if(type == 1) return -1.5 * (120.0) * std::exp(-120.0 * c)
        + (0.0351 / 0.083) * (1.0 / (std::cosh((c - 0.286) / 0.083) * std::cosh((c - 0.286) / 0.083)))
        - (0.0045 / 0.119) * (1.0 / (std::cosh((c - 0.849) / 0.119) * std::cosh((c - 0.849) / 0.119)))
        - (0.035 / 0.05) * (1.0 / (std::cosh((c - 0.9233) / 0.05) * std::cosh((c - 0.9233) / 0.05)))
        - (0.0147 / 0.034) * (1.0 / (std::cosh((c - 0.5) / 0.034) * std::cosh((c - 0.5) / 0.034)))
        - (0.102 / 0.142) * (1.0 / (std::cosh((c - 0.194) / 0.142) * std::cosh((c - 0.194) / 0.142)))
        - (0.022 / 0.0164) * (1.0 / (std::cosh((c - 0.9) / 0.0164) * std::cosh((c - 0.9) / 0.0164)))
        - (0.011 / 0.0226) * (1.0 / (std::cosh((c - 0.124) / 0.0226) * std::cosh((c - 0.124) / 0.0226)))
        + (0.0155 / 0.029) * (1.0 / (std::cosh((c - 0.105) / 0.029) * std::cosh((c - 0.105) / 0.029)));
    else return 0.07645 * (-54.4806) * (1.0 / (std::cosh(30.834 - 54.4806 * c) * std::cosh(30.834 - 54.4806 * c)))
        + 2.1581 * (-50.294) * (1.0 / (std::cosh(52.294 - 50.294 * c) * std::cosh(52.294 - 50.294 * c)))
        - 0.14169 * (-19.8543) * (1.0 / (std::cosh(11.0923 - 19.8543 * c) * std::cosh(11.0923 - 19.8543 * c)))
        + 0.2051 * (-5.4888) * (1.0 / (std::cosh(1.4684 - 5.4888 * c) * std::cosh(1.4684 - 5.4888 * c)))
        - 0.2531 / (0.1316) * (1.0 / (std::cosh((-c + 0.56468) / 0.1316) * std::cosh((-c + 0.56468) / 0.1316)))
        - 0.02167 / (0.006) * (1.0 / (std::cosh((c - 0.525) / 0.006) * std::cosh((c - 0.525) / 0.006)));
}

double dc_dx(double c) {
    return std::sqrt(c * (1.0 - c)) * 2 * constant::k;
}

double dc_dxdc(double c) {
    return ((1.0 - 2.0 * c) / (2 * std::sqrt(c * (1.0 - c)))) * 2 * constant::k;
}

double j_border(double c, int type) {
    return -dc_dx(c) * std::sinh((constant::delta_u - ueq(c, type)) * 0.5);
}

double d_j_border(double c, int type) {
    return -dc_dxdc(c) * ueq(c, type) + 0.5 * dc_dx(c) * d_ueq(c, type) * std::cosh((constant::delta_u - ueq(c, type)) * 0.5);
}

void particle_solver::calc(Eigen::Ref<MatrixXd> u) {
    long element_cnt = this->point_coord.size() - 1;
    result = VectorXd::Zero(element_cnt + 1);

    int iter_time = 0;
    double ratio = 999999.0;
    double first_norm;
    double eff = M_PI * 4.0;
    MatrixXd du = MatrixXd::Zero(element_cnt + 1, 1);

    while(iter_time < iter && ratio > tolerance) {
        this->gen.generate(u, du);
            
        gen.generated_K(element_cnt, element_cnt) -= d_j_border(u(element_cnt, 0), type) * eff;
        gen.generated_res(element_cnt) -= j_border(u(element_cnt, 0), type) * eff;

        MatrixXd delta = -gen.generated_K.lu().solve(gen.generated_res);
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
