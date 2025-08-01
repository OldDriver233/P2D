#include "stiffness_stress.h"


void stiffness_stress::generate(std::vector<Eigen::Triplet<double> > &t, std::vector<Eigen::Triplet<double>> &l) {
    const int dim = get_dim(mesh.p_type), n = get_nodes(mesh.p_type), dim_voigt = (dim * (dim + 1)) / 2;
    std::size_t node_cnt = mesh.node_count;

    MatrixXd xs;
    MatrixXd w;

    if (dim == 1) {
        xs = get_integration_point<1, 2>();
        w = get_integration_weight<1, 2>();
    } else if (dim == 2) {
        if (n == 3) {
            xs = get_integration_point<2, 3>();
            w = get_integration_weight<2, 3>();
        } else {
            xs = get_integration_point<2, 4>();
            w = get_integration_weight<2, 4>();
        }
    }

    for (int e = 0; e < mesh.elem_count; e++) {
        MatrixXd B(dim_voigt, dim * n), K(dim * n, dim * n), C(dim_voigt, dim_voigt);
        K = MatrixXd::Zero(dim * n, dim * n);
        double E = 70e9, nu = .26;
        double E_ref = 70e9;
        C << 1, nu, 0,
             nu, 1, 0,
             0, 0, (1 - nu) / 2;
        C = C * E / (1 - nu * nu);

        for (int j = 0; j < n; j++) {
            const MatrixXd &Ns = shapes.cached_matrix_N[e * n + j];
            const MatrixXd &dNs = shapes.cached_matrix_dN[e * n + j];
            double det = shapes.cached_det_J[e * n + j];

            if (dim == 2) {
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < n; l++) {
                        B(k, l * dim + k) = dNs(l, k);
                    }
                }
                for (int k = 0; k < dim * n; k++) {
                    B(2, k) = dNs(k / 2, (k + 1) % 2);
                }
            }

            K += B.transpose() * C * B * w(j) * det / E_ref;
        }

        //std::cout<<K<<std::endl<<std::endl;

        for (int j = 0; j < n * dim; j++) {
            std::size_t id_l = mesh.elements[e * n + j / dim];
            std::size_t dof_l = 5 + j % dim;
            for (int k = 0; k < n * dim; k++) {
                std::size_t id_r = mesh.elements[e * n + k / dim];
                std::size_t dof_r = 5 + k % dim;
                if (!mesh.anode_cc_wall_nodes.contains(id_l) && !mesh.anode_cc_wall_nodes.contains(id_r)) {
                    t.emplace_back(dof.get_dof(id_l, dof_l), dof.get_dof(id_r, dof_r), K(j, k));
                    l.emplace_back(id_l * dim + dof_l - 5, id_r * dim + dof_r - 5, K(j, k));
                }
            }
            if (mesh.anode_cc_wall_nodes.contains(id_l)) {
                t.emplace_back(dof.get_dof(id_l, dof_l), dof.get_dof(id_l, dof_l), 1);
                l.emplace_back(id_l * dim + dof_l - 5, id_l * dim + dof_l - 5, 1);
            }
        }
    }
}

void stiffness_stress::generate_residue(const Eigen::Ref<MatrixXd> &u, const std::vector<double> avg_c_an,
    const std::vector<double> avg_c_ca, Eigen::Ref<VectorXd> res, const Eigen::Ref<Eigen::SparseMatrix<double>>& K) {
    const int dim = get_dim(mesh.p_type), n = get_nodes(mesh.p_type), dim_voigt = (dim * (dim + 1)) / 2;
    std::size_t node_cnt = mesh.node_count;
    double E_ref = 70e9;

    MatrixXd xs;
    MatrixXd w;

    if (dim == 1) {
        xs = get_integration_point<1, 2>();
        w = get_integration_weight<1, 2>();
    } else if (dim == 2) {
        if (n == 3) {
            xs = get_integration_point<2, 3>();
            w = get_integration_weight<2, 3>();
        } else {
            xs = get_integration_point<2, 4>();
            w = get_integration_weight<2, 4>();
        }
    }

    VectorXd load = VectorXd::Zero(dim * node_cnt);
    VectorXd local_u(dim * node_cnt);

    for (int x = 0; x < node_cnt; x++) {
        for (int j = 0; j < dim; j++) {
            local_u(x * dim + j) = u(dof.get_dof(x, j + 5), 0);
        }
    }


    for (int e = 0; e < mesh.elem_count; e++) {
        VectorXd e_load = VectorXd::Zero(n * dim);
        for (int j = 0; j < n; j++) {
            const MatrixXd &N = shapes.cached_matrix_N_mult[e * n + j];
            const double det = shapes.cached_det_J[e * n + j];
            VectorXd F(2);
            F<<1000, 0;
            e_load += N.transpose() * F * w(j) * det / E_ref;
        }


        for (int j = 0; j < n * dim; j++) {
            std::size_t id_l = mesh.elements[e * n + j / dim];
            std::size_t dof_l = 5 + j % dim;
            if (!mesh.anode_cc_wall_nodes.contains(id_l)) {
                load(id_l * dim + dof_l - 5) += e_load(j);
            }
        }
    }

    VectorXd local_res = K * local_u - load;
    for (int x = 0; x < node_cnt; x++) {
        for (int j = 0; j < dim; j++) {
            res(dof.get_dof(x, j + 5)) = local_res(x * dim + j);
        }
    }
}
