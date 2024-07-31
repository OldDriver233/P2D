#include "stiffness_generator.h"
#include "../integration/integration_shapes.h"
#include "../shaping/primitives.h"
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>
#include <cmath>


stiffness_generator::stiffness_generator(VectorXd points, Primitive p, double dt, double d_ref, double d, double r) : points(
        points), dt(dt), d_ref(d_ref), d(d), r(r) {
    switch (p) {
        case Primitive::Line2:
            dim = 1;
            n = 2;
            break;
        case Primitive::Line3:
            dim = 1;
            n = 3;
            break;
    }

    // Next, generate N and dN
    long elem_cnt = points.size() - 1;
    MatrixXd x = get_integration_point(dim, n);
    for (auto i = 0; i < elem_cnt; i++) {
        MatrixXd N(n, 1), dN(n, 1);
        MatrixXd coords(1, n);
        coords << points(i), points(i + 1);
        for (auto j = 0; j < n; j++) {
            N = get_shape_func_at(dim, n, x(j));
            auto dNds = get_shape_deriv_at(dim, n, x(j));

            // NOTE: 1D situation, J is 1x1
            auto J = coords * dNds;
            double det_J = J.determinant();
            dN = dNds / det_J;

            this->cached_matrix_N.push_back(N);
            this->cached_matrix_dN.push_back(dN);
            this->cached_det_J.push_back(det_J);
        }
    }
}

void stiffness_generator::generate(Eigen::Ref<MatrixXd> u, Eigen::Ref<MatrixXd> du) {
    long elem_cnt = this->points.size() - 1;
    this->generated_K = MatrixXd::Zero(elem_cnt + 1, elem_cnt + 1);
    this->generated_res = VectorXd::Zero(elem_cnt + 1);
    MatrixXd xs = get_integration_point(dim, n);
    MatrixXd w = get_integration_weight(dim, n);
    double d_s = d / d_ref;
    double eff_1 = (r * r / d_ref) / dt;
    double eff_2 = d_s;

    for (int i = 0; i < elem_cnt; i++) {
        MatrixXd e_u = u({i, i + 1}, 0);
        MatrixXd e_du = du({i, i + 1}, 0);
        MatrixXd e_k = MatrixXd::Zero(n, n);
        MatrixXd e_r = MatrixXd::Zero(n, 1);
        for (int j = 0; j < n; j++) {
            const MatrixXd &N = cached_matrix_N[i * n + j];
            const MatrixXd &dN = cached_matrix_dN[i * n + j];
            MatrixXd N_T = N.transpose();
            MatrixXd dN_T = dN.transpose();
            double det = cached_det_J[i * n + j];
            double s = xs(j);
            double lower = this->points[i];
            double upper = this->points[i + 1];
            double x = lower + (s + 1) * (upper - lower) / 2;

            e_r += N * N_T * e_du * 4 * M_PI * x * x * w(j) * det * eff_1;
            e_r += dN * dN_T * e_u * 4 * M_PI * x * x * w(j) * det * eff_2;
            e_k += N * N_T * 4 * M_PI * x * x * w(j) * det * eff_1;
            e_k += dN * dN_T * 4 * M_PI * x * x * w(j) * det * eff_2;

        }

        for (int j = 0; j < n; j++) {
            for (int l = 0; l < n; l++) {
                this->generated_K(i + j, i + l) += e_k(j, l);
            }
            this->generated_res(i + j) += e_r(j);
        }
    }
}