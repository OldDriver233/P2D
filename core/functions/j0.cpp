#include "functions.h"
#include <cmath>

double j0(double c_e, double c_a, int type) {
    const double k = (type == 1 ? 5.031e-11 : 2.334e-11);
    const double c_max = (type == 1 ? 30555 : 51554);
    return k * std::sqrt(c_e * (c_max - c_a) * c_a);
}

double d_j0_e(double c_e, double c_a, int type) {
    const double k = (type == 1 ? 5.031e-11 : 2.334e-11);
    const double c_max = (type == 1 ? 30555 : 51554);
    return k * std::sqrt(c_a * (c_max - c_a)) / std::sqrt(c_e) * 0.5;
}

double d_j0_a(double c_e, double c_a, int type) {
    const double k = (type == 1 ? 5.031e-11 : 2.334e-11);
    const double c_max = (type == 1 ? 30555 : 51554);
    return -k * std::sqrt(c_e) * std::sqrt(c_a) * 0.5 / std::sqrt(c_max - c_a)
            + k * std::sqrt(c_max - c_a) * std::sqrt(c_e) * 0.5 / std::sqrt(c_a);
}