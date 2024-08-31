#ifndef FEM_STIFFNESS_BASE_H
#define FEM_STIFFNESS_BASE_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <vector>
#include "../../integration/integration_shapes.h"
#include "../../shaping/primitives.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_base{
public:
    std::vector<MatrixXd> cached_matrix_N, cached_matrix_dN;
    std::vector<double> cached_det_J;
    MatrixXd points;
    int surface_an_sep, surface_ca_sep;

    stiffness_base() {}
    ~stiffness_base() {
        this->cached_matrix_N.clear();
        this->cached_matrix_dN.clear();
        this->cached_det_J.clear();
    }
    stiffness_base(VectorXd points): points(points) {
        int dim = 1, n = 2;
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

    virtual void generate(Eigen::Ref<MatrixXd>, Eigen::Ref<MatrixXd>, Eigen::Ref<MatrixXd>, Eigen::Ref<VectorXd>) = 0;
};

#endif //FEM_STIFFNESS_BASE_H