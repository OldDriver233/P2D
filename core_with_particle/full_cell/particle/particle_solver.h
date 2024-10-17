#ifndef FEM_PARTICLE_SOLVER_H
#define FEM_PARTICLE_SOLVER_H
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <iostream>
#include <ostream>
#include "../../constants/constant.h"
#include "../../integration/integration_shapes.h"
#include "../../shaping/primitives.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

class particle_solver {
public:
    VectorXd point_coord;
    Eigen::SparseMatrix<double> assembled_A;
    Eigen::SparseMatrix<double> assembled_B;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_A;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_B;
    MatrixXd last_cs;
    VectorXd j_coeff;
    const double D_s;
    const double D_sref;
    const double c_max;

    particle_solver(const VectorXd& coord, const double D_s, const double c_max): point_coord(coord), D_s(D_s), D_sref(D_s), c_max(c_max) {
        const int dim = 1, n = 2;
        const double R_s = constant::r_p;
        int pt_size = coord.size();
        int elem_size = pt_size - 1;

        assembled_A = Eigen::SparseMatrix<double>(pt_size, pt_size);
        assembled_B = Eigen::SparseMatrix<double>(pt_size, pt_size);
        std::vector<Eigen::Triplet<double>> coeff_A;
        std::vector<Eigen::Triplet<double>> coeff_B;
        j_coeff = VectorXd::Zero(pt_size);
        j_coeff(pt_size - 1) = R_s / D_sref * 4 * M_PI * constant::j_ref / c_max;

        MatrixXd xs = get_integration_point<dim, n>();
        MatrixXd w = get_integration_weight<dim, n>();

        const double eff_1 = R_s * R_s / D_sref / constant::dt * 4 * M_PI;
        const double eff_2 = D_s / D_sref * 4 * M_PI;

        for(int i = 0; i < elem_size; i++) {
            MatrixXd coords(1, n);
            coords << coord(i), coord(i + 1);
            MatrixXd e_a = MatrixXd::Zero(n, n);
            MatrixXd e_b = MatrixXd::Zero(n, n);
            for(int j = 0; j < n; j++) {
                MatrixXd N = get_shape_func_at<dim, n>(xs(j));
                VectorXd dNds = get_shape_deriv_at<dim, n>(xs(j));
                VectorXd J = coords * dNds;
                double det_J = J.determinant();
                MatrixXd dN = dNds / det_J;
                MatrixXd N_T = N.transpose();
                MatrixXd dN_T = dN.transpose();


                double s = xs(j);
                double lower = coords(0, 0);
                double upper = coords(0, 1);
                double x = lower + (s + 1) * (upper - lower) / 2;

                e_a += N * N_T * x * x * w(j) * det_J * eff_1 + dN * dN_T * x * x * w(j) * det_J * eff_2;
                e_b += N * N_T * x * x * w(j) * det_J * eff_1;
            }

            for(int j = 0; j < n; j++) {
                for(int l = 0; l < n; l++) {
                    coeff_A.emplace_back(i + j, i + l, e_a(j, l));
                    coeff_B.emplace_back(i + j, i + l ,e_b(j, l));
                }
            }
        }

        assembled_A.setFromTriplets(coeff_A.begin(), coeff_A.end());
        assembled_B.setFromTriplets(coeff_B.begin(), coeff_B.end());
        solver_A.compute(assembled_A);
        solver_B.compute(assembled_B);
        j_coeff = solver_A.solve(j_coeff);
    }

    void pre_calc(const Eigen::Ref<MatrixXd> &c_s);
    void calc(Eigen::Ref<MatrixXd> c_s, const Eigen::Ref<MatrixXd> &u, int pt_size, int an, int ca, int type);
};

#endif //FEM_PARTICLE_SOLVER_H