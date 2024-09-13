#include "functions.h"
#include <cmath>

double j0(double c_e, double c_a, int type) {
    const double k = (type == 1 ? 1e-10 : 3e-11);
    return k * std::sqrt(c_e * (1 - c_a) * c_a);
}

double d_j0_e(double c_e, double c_a, int type) {
    const double k = (type == 1 ? 1e-10 : 3e-11);
    return k * std::sqrt(c_a * (1 - c_a)) / std::sqrt(c_e) * 0.5;
}

double d_j0_a(double c_e, double c_a, int type) {
    const double k = (type == 1 ? 1e-10 : 3e-11);
    return -k * std::sqrt(c_e) * std::sqrt(c_a) * 0.5 / std::sqrt(1 - c_a)
            + k * std::sqrt(1 - c_a) * std::sqrt(c_e) * 0.5 / std::sqrt(c_a);
}