#include "particle_solver.h"
#include "stiffness_generator.h"
#include <iostream>

void particle_solver::calc(Eigen::Ref<MatrixXd> u) {
    long element_cnt = this->point_coord.size() - 1;
    result = VectorXd::Zero(element_cnt + 1);

    int iter_time = 0;
    double ratio = 999999.0;
    double first_norm;
    MatrixXd du = MatrixXd::Zero(element_cnt + 1, 1);

    while(iter_time < iter && ratio > tolerance) {
        this->gen.generate(u, du);
        gen.generated_res(element_cnt) -= 1e-9 * M_PI * 4.0 * r / d_ref;

        MatrixXd delta = -gen.generated_K.llt().solve(gen.generated_res);
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
