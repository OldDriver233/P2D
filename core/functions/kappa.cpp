#include "functions.h"

double kappa(double c) {
    const double a_0 = 0.0911;
    const double a1 = 1.9101e-3;
    const double a2 = -1.052e-6;
    const double a3 = 0.1554e-9;
    return a_0 + a1 * c + a2 * c * c + a3 * c * c * c;
}

double d_kappa(double c) {
    const double a_0 = 0.0911;
    const double a1 = 1.9101e-3;
    const double a2 = -1.052e-6;
    const double a3 = 0.1554e-9;
    return a1 + 2 * a2 * c + 3 * a3 * c * c;
}