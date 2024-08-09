#include "particle_solver.h"
#include "stiffness_generator.h"
#include <cmath>
#include <eigen3/Eigen/src/Core/ArithmeticSequence.h>
#include <iostream>

const double delta_u = 0.04;
const double max_c = 28740.0;

double clamp(double x, double lower, double upper) {
    return x < lower ? lower : (x > upper ? upper : x);
}

double ueq(double c) {
    return 0.194 + 1.5 * std::exp(-120.0 * c)
        + 0.0351 * std::tanh((c - 0.286) / 0.083)
        - 0.0045 * std::tanh((c - 0.849) / 0.119)
        - 0.035 * std::tanh((c - 0.9233) / 0.05)
        - 0.0147 * std::tanh((c - 0.5) / 0.034)
        - 0.102 * std::tanh((c - 0.194) / 0.142)
        - 0.022 * std::tanh((c - 0.9) / 0.0164)
        - 0.011 * std::tanh((c - 0.124) / 0.0226)
        + 0.0155 * std::tanh((c - 0.105) / 0.029);
}

double d_ueq(double c) {
    return -1.5 * (120.0 / max_c) * std::exp(-120.0 * c)
        + (0.0351 / (0.083 * max_c)) * (1.0 / (std::cosh((c - 0.286) / 0.083) * std::cosh((c - 0.286) / 0.083)))
        - (0.0045 / (0.119 * max_c)) * (1.0 / (std::cosh((c - 0.849) / 0.119) * std::cosh((c - 0.849) / 0.119)))
        - (0.035 / (0.05 * max_c)) * (1.0 / (std::cosh((c - 0.9233) / 0.05) * std::cosh((c - 0.9233) / 0.05)))
        - (0.0147 / (0.034 * max_c)) * (1.0 / (std::cosh((c - 0.5) / 0.034) * std::cosh((c - 0.5) / 0.034)))
        - (0.102 / (0.142 * max_c)) * (1.0 / (std::cosh((c - 0.194) / 0.142) * std::cosh((c - 0.194) / 0.142)))
        - (0.022 / (0.0164 * max_c)) * (1.0 / (std::cosh((c - 0.9) / 0.0164) * std::cosh((c - 0.9) / 0.0164)))
        - (0.011 / (0.0226 * max_c)) * (1.0 / (std::cosh((c - 0.124) / 0.0226) * std::cosh((c - 0.124) / 0.0226)))
        + (0.0155 / (0.029 * max_c)) * (1.0 / (std::cosh((c - 0.105) / 0.029) * std::cosh((c - 0.105) / 0.029)));
}

double dc_dx(double c) {
    return std::sqrt(c * (1.0 - c)) * 0.3478;
}

double dc_dxdc(double c) {
    return ((1.0 - 2.0 * c) / (2 * std::sqrt(c * (1.0 - c)))) * 0.3478;
}

double j_border(double c) {
    return dc_dx(c) * std::sinh(delta_u - ueq(c));
}

double d_j_border(double c) {
    return dc_dxdc(c) * ueq(c) - dc_dx(c) * d_ueq(c) * std::cosh(delta_u - ueq(c));
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
            
        gen.generated_K(element_cnt, element_cnt) -= d_j_border(u(element_cnt, 0)) * eff;
        gen.generated_res(element_cnt) -= j_border(u(element_cnt, 0)) * eff;

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
