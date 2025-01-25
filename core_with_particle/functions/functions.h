#ifndef FEM_FUNCTIONS_H
#define FEM_FUNCTIONS_H
#include "../constants/constant.h"
#include <cmath>
#include <numbers>
#include <eigen3/Eigen/Dense>

using Eigen::VectorXd;

inline double bv(double eta) {
    const double F = 96485.3329;
    const double R = 8.3145;
    const double T = 298.15;

    return 2 * std::sinh(0.5 * F * eta / (R * T));
}

inline double d_bv(double eta) {
    const double F = 96485.3329;
    const double R = 8.3145;
    const double T = 298.15;

    return 2 * 0.5 * F / (R * T) * std::cosh(0.5 * F * eta / (R * T));
}

template<int type>
inline double j0(double c_e, double c_a, double t = constant::t_ref) {
    const double k = (type == 1 ? constant::k_an : constant::k_ca);
    const double E = (type == 1 ? constant::exchange_energy_an : constant::exchange_energy_ca);
    const double k_ref = k * exp((1 / t - 1 / constant::t_ref) * -E / constant::R);
    return k_ref * std::sqrt(c_e * (1 - c_a) * c_a);
}

template<int type>
inline double d_j0_e(double c_e, double c_a, double t = constant::t_ref) {
    const double k = (type == 1 ? constant::k_an : constant::k_ca);
    const double E = (type == 1 ? constant::exchange_energy_an : constant::exchange_energy_ca);
    const double k_ref = k * exp((1 / t - 1 / constant::t_ref) * -E / constant::R);
    return k_ref * std::sqrt(c_a * (1 - c_a)) / std::sqrt(c_e) * 0.5;
}

template<int type>
inline double d_j0_a(double c_e, double c_a, double t = constant::t_ref) {
    const double k = (type == 1 ? constant::k_an : constant::k_ca);
    const double E = (type == 1 ? constant::exchange_energy_an : constant::exchange_energy_ca);
    const double k_ref = k * exp((1 / t - 1 / constant::t_ref) * -E / constant::R);
    return -k_ref * std::sqrt(c_e) * std::sqrt(c_a) * 0.5 / std::sqrt(1 - c_a) +
           k_ref * std::sqrt(1 - c_a) * std::sqrt(c_e) * 0.5 / std::sqrt(c_a);
}

template<int type>
inline double d_j0_t(double c_e, double c_a, double t = constant::t_ref) {
    const double k = (type == 1 ? constant::k_an : constant::k_ca);
    const double E = (type == 1 ? constant::exchange_energy_an : constant::exchange_energy_ca);
    const double k_ref = k * exp((1 / t - 1 / constant::t_ref) * -E / constant::R);
    return E / constant::R * k_ref / (t * t);
}

inline double kappa(double c) {
    const double a_0 = 0.0911;
    const double a1 = 1.9101e-3;
    const double a2 = -1.052e-6;
    const double a3 = 0.1554e-9;
    return a_0 + a1 * c + a2 * c * c + a3 * c * c * c;
}

inline double d_kappa(double c) {
    const double a_0 = 0.0911;
    const double a1 = 1.9101e-3;
    const double a2 = -1.052e-6;
    const double a3 = 0.1554e-9;
    return a1 + 2 * a2 * c + 3 * a3 * c * c;
}

template<int type>
inline double uoc(double c) {
    if (type == 1)
        return 0.194 + 1.5 * std::exp(-120.0 * c) +
               0.0351 * std::tanh((c - 0.286) / 0.083) -
               0.0045 * std::tanh((c - 0.849) / 0.119) -
               0.035 * std::tanh((c - 0.9233) / 0.05) -
               0.0147 * std::tanh((c - 0.5) / 0.034) -
               0.102 * std::tanh((c - 0.194) / 0.142) -
               0.022 * std::tanh((c - 0.9) / 0.0164) -
               0.011 * std::tanh((c - 0.124) / 0.0226) +
               0.0155 * std::tanh((c - 0.105) / 0.029);
    else {
        c *= 1.062;
        return 2.16216 + 0.07645 * std::tanh(30.834 - 54.4806 * c) +
               2.1581 * std::tanh(52.294 - 50.294 * c) -
               0.14169 * std::tanh(11.0923 - 19.8543 * c) +
               0.2051 * std::tanh(1.4684 - 5.4888 * c) +
               0.2531 * std::tanh((-c + 0.56478) / 0.1316) -
               0.02167 * std::tanh((c - 0.525) / 0.006);
    }
}

template<int type>
inline double d_uoc(double c) {
    if (type == 1) {
        return -1.5 * (120.0) * std::exp(-120.0 * c) +
               (0.0351 / 0.083) * (1.0 / (std::cosh((c - 0.286) / 0.083) *
                                          std::cosh((c - 0.286) / 0.083))) -
               (0.0045 / 0.119) * (1.0 / (std::cosh((c - 0.849) / 0.119) *
                                          std::cosh((c - 0.849) / 0.119))) -
               (0.035 / 0.05) * (1.0 / (std::cosh((c - 0.9233) / 0.05) *
                                        std::cosh((c - 0.9233) / 0.05))) -
               (0.0147 / 0.034) * (1.0 / (std::cosh((c - 0.5) / 0.034) *
                                          std::cosh((c - 0.5) / 0.034))) -
               (0.102 / 0.142) * (1.0 / (std::cosh((c - 0.194) / 0.142) *
                                         std::cosh((c - 0.194) / 0.142))) -
               (0.022 / 0.0164) * (1.0 / (std::cosh((c - 0.9) / 0.0164) *
                                          std::cosh((c - 0.9) / 0.0164))) -
               (0.011 / 0.0226) * (1.0 / (std::cosh((c - 0.124) / 0.0226) *
                                          std::cosh((c - 0.124) / 0.0226))) +
               (0.0155 / 0.029) * (1.0 / (std::cosh((c - 0.105) / 0.029) *
                                          std::cosh((c - 0.105) / 0.029)));
    } else {
        c *= 1.062;
        return 0.07645 * (-54.4806) *
               (1.0 / (std::cosh(30.834 - 54.4806 * c) *
                       std::cosh(30.834 - 54.4806 * c))) +
               2.1581 * (-50.294) *
               (1.0 / (std::cosh(52.294 - 50.294 * c) *
                       std::cosh(52.294 - 50.294 * c))) -
               0.14169 * (-19.8543) *
               (1.0 / (std::cosh(11.0923 - 19.8543 * c) *
                       std::cosh(11.0923 - 19.8543 * c))) +
               0.2051 * (-5.4888) *
               (1.0 / (std::cosh(1.4684 - 5.4888 * c) *
                       std::cosh(1.4684 - 5.4888 * c))) -
               0.2531 / (0.1316) *
               (1.0 / (std::cosh((-c + 0.56468) / 0.1316) *
                       std::cosh((-c + 0.56468) / 0.1316))) -
               0.02167 / (0.006) *
               (1.0 / (std::cosh((c - 0.525) / 0.006) *
                       std::cosh((c - 0.525) / 0.006)));
    }
}
#endif // FEM_FUNCTIONS_H
