#ifndef PRE_CALC_H
#define PRE_CALC_H
#include <vector>
#include <Eigen/Dense>

#include "../../integration/integration_shapes.h"
#include "../../mesh/mesh_reader.h"
#include "../../mesh/dof_assigner.h"
#include "../../shaping/primitives.h"

using Eigen::MatrixXd;

struct pre_calc_shapes {
    std::vector<MatrixXd> cached_matrix_N, cached_matrix_dN;
    std::vector<MatrixXd> cached_matrix_NNT, cached_matrix_dNdNT;
    std::vector<MatrixXd> cached_matrix_NdNT;
    std::vector<double> cached_det_J;
    const mesh_reader& mesh;
    const dof_assigner& dof;

    pre_calc_shapes(const mesh_reader& mesh, const dof_assigner& dof): mesh(mesh), dof(dof) {
        const int dim = get_dim(mesh.p_type);
        const int n = get_nodes(mesh.p_type);
        MatrixXd x;
        if (dim == 1) {
            x = get_integration_point<1, 2>();
        } else if (dim == 2) {
            if (n == 3) x = get_integration_point<2, 3>();
            else x = get_integration_point<2, 4>();
        }
        for (auto i = 0; i < mesh.elem_count; i++) {
            MatrixXd N(n, 1), dN(n, dim), dNdu;
            MatrixXd coord(dim, n);
            for (auto j = 0; j < n; j++) {
                coord.col(j) = mesh.coord.col(mesh.elements[i * n + j]);
            }
            for (auto j = 0; j < n; j++) {
                if (dim == 1) {
                    N = get_shape_func_at<1, 2>(x.row(j).transpose());
                    dNdu = get_shape_deriv_at<1, 2>(x.row(j).transpose());
                } else if (dim == 2) {
                    if (n == 3) {
                        N = get_shape_func_at<2, 3>(x.row(j).transpose());
                        dNdu = get_shape_deriv_at<2, 3>(x.row(j).transpose());
                    } else {
                        N = get_shape_func_at<2, 4>(x.row(j).transpose());
                        dNdu = get_shape_deriv_at<2, 4>(x.row(j).transpose());
                    }
                }

                MatrixXd J = coord * dNdu;
                dN = dNdu * J.inverse();

                this->cached_matrix_N.push_back(N);
                this->cached_matrix_dN.push_back(dN);
                this->cached_matrix_NdNT.push_back(N * dN.transpose());
                this->cached_matrix_NNT.push_back(N * N.transpose());
                this->cached_matrix_dNdNT.push_back(dN * dN.transpose());
                this->cached_det_J.push_back(J.determinant());
            }
        }
    }
};

#endif //PRE_CALC_H
