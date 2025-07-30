#ifndef FEM_PARTICLE_SOLVER_H
#define FEM_PARTICLE_SOLVER_H
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <iostream>
#include <ostream>
#include "../../constants/constant.h"
#include "../../integration/integration_shapes.h"
#include "../../io/settings/settings.h"
#include "../../shaping/primitives.h"
#include "../../mesh/mesh_reader.h"
#include "../../mesh/dof_assigner.h"

class StepControl;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class particle_solver {
public:
    VectorXd point_coord;
    Eigen::SparseMatrix<double> assembled_A_1;
    Eigen::SparseMatrix<double> assembled_A_2;
    Eigen::SparseMatrix<double> assembled_B;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> constant_solver;
    MatrixXd last_cs;
    VectorXd pre_j_coeff;
    VectorXd j_coeff;
    const double D_s;
    const double D_sref;
    const double c_max;
    StepControl* step_control;

    particle_solver(const VectorXd& coord, const double D_s, const double c_max, StepControl* st): point_coord(coord), D_s(D_s), D_sref(D_s), c_max(c_max), step_control(st) {
        const int dim = 1, n = 2;
        const double R_s = constant::r_p;
        int pt_size = coord.size();
        int elem_size = pt_size - 1;

        assembled_A_1 = Eigen::SparseMatrix<double>(pt_size, pt_size);
        assembled_A_2 = Eigen::SparseMatrix<double>(pt_size, pt_size);
        assembled_B = Eigen::SparseMatrix<double>(pt_size, pt_size);
        std::vector<Eigen::Triplet<double>> coeff_A_1;
        std::vector<Eigen::Triplet<double>> coeff_A_2;
        std::vector<Eigen::Triplet<double>> coeff_B;
        pre_j_coeff = VectorXd::Zero(pt_size);
        if (!settings::stress_analysis) pre_j_coeff(pt_size - 1) = R_s / D_sref * 4 * M_PI * constant::j_ref / c_max;
        else pre_j_coeff(pt_size - 1) = R_s / D_sref * constant::j_ref / c_max;

        MatrixXd xs = get_integration_point<dim, n>();
        MatrixXd w = get_integration_weight<dim, n>();

        const double eff_1 = R_s * R_s / D_sref * 4 * M_PI;
        const double eff_2 = D_s / D_sref * 4 * M_PI;

        for(int i = 0; i < elem_size; i++) {
            MatrixXd coords(1, n);
            coords << coord(i), coord(i + 1);
            MatrixXd e_a_1 = MatrixXd::Zero(n, n);
            MatrixXd e_a_2 = MatrixXd::Zero(n, n);
            MatrixXd e_b = MatrixXd::Zero(n, n);
            for(int j = 0; j < n; j++) {
                MatrixXd N = get_shape_func_at<dim, n>(xs.row(j).transpose());
                VectorXd dNds = get_shape_deriv_at<dim, n>(xs.row(j).transpose());
                VectorXd J = coords * dNds;
                double det_J = J.determinant();
                MatrixXd dN = dNds * J.inverse();
                MatrixXd N_T = N.transpose();
                MatrixXd dN_T = dN.transpose();


                double s = xs(j);
                double lower = coords(0, 0);
                double upper = coords(0, 1);
                double x = lower + (s + 1) * (upper - lower) / 2;

                e_a_1 += dN * dN_T * x * x * w(j) * det_J * eff_2;
                e_a_2 += N * N_T * x * x * w(j) * det_J * eff_1;
                e_b += N * N_T * x * x * w(j) * det_J * eff_1;
            }

            for(int j = 0; j < n; j++) {
                for(int l = 0; l < n; l++) {
                    coeff_A_1.emplace_back(i + j, i + l, e_a_1(j, l));
                    coeff_A_2.emplace_back(i + j, i + l, e_a_2(j, l));
                    coeff_B.emplace_back(i + j, i + l ,e_b(j, l));
                }
            }
        }

        assembled_A_1.setFromTriplets(coeff_A_1.begin(), coeff_A_1.end());
        assembled_A_2.setFromTriplets(coeff_A_2.begin(), coeff_A_2.end());
        assembled_B.setFromTriplets(coeff_B.begin(), coeff_B.end());

        constant_solver.compute(assembled_A_1 + assembled_A_2 / constant::dt);
        j_coeff = constant_solver.solve(pre_j_coeff);
    }

    void pre_calc(const Eigen::Ref<MatrixXd> &c_s);
    double calc(Eigen::Ref<MatrixXd> c_s, const Eigen::Ref<MatrixXd> &u, int type, const mesh_reader& mesh, const dof_assigner& dof);
    double calc_stress(Eigen::Ref<MatrixXd> c_s, const Eigen::Ref<MatrixXd> &u, int type, const mesh_reader& mesh, const dof_assigner& dof);
};

#endif //FEM_PARTICLE_SOLVER_H