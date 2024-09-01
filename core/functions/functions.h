#ifndef FEM_FUNCTIONS_H
#define FEM_FUNCTIONS_H
#include <cmath>
#include <numbers>

double uoc(double c, int type);
double d_uoc(double c, int type);
double j0(double c_e, double c_a, int type);
double d_j0_e(double c_e, double c_a, int type);
double d_j0_a(double c_e, double c_a, int type);
double bv(double eta);
double d_bv(double eta);

#endif //FEM_FUNCTIONS_H