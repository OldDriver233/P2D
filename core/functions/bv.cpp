#include "functions.h"

double bv(double eta) {
    const double F = 96485.3329;
    const double R = 8.3145;
    const double T = 298.15;

    return 2 * std::sinh(0.5 * F * eta / (R * T));
}

double d_bv(double eta) {
    const double F = 96485.3329;
    const double R = 8.3145;
    const double T = 298.15;

    return 2 * 0.5 * F / (R * T) * std::cosh(0.5 * F * eta / (R * T));
}