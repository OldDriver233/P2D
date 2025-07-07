#include "stiffness_stress.h"

#include "../../constants/constant.h"

std::tuple<double, double> stiffness_stress::get_material_property(std::size_t element) const {
    double E, nu;
    if (mesh.anode_element_set.contains(element)) {
        E = constant::anode.E;
        nu = constant::anode.nu;
    } else if (mesh.cathode_element_set.contains(element)) {
        E = constant::cathode.E;
        nu = constant::cathode.nu;
    } else if (mesh.separator_element_set.contains(element)) {
        E = constant::separator.E;
        nu = constant::separator.nu;
    } else if (mesh.anode_cc_element_set.contains(element)) {
        E = constant::anode_cc.E;
        nu = constant::anode_cc.nu;
    } else if (mesh.cathode_cc_element_set.contains(element)) {
        E = constant::cathode_cc.E;
        nu = constant::cathode_cc.nu;
    }
    return std::make_tuple(E, nu);
}


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
        MatrixXd K(dim * n, dim * n), C(dim_voigt, dim_voigt);
        K = MatrixXd::Zero(dim * n, dim * n);
        double E_ref = 70e9;
        auto [E, nu] = this->get_material_property(e);
        C << 1 - nu, nu, 0,
             nu, 1 - nu, 0,
             0, 0, (1 - 2 * nu) / 2;
        C = C * E / (1 + nu) / (1 - 2 * nu);

        for (int j = 0; j < n; j++) {
            const MatrixXd &Ns = shapes.cached_matrix_N[e * n + j];
            const MatrixXd &dNs = shapes.cached_matrix_dN[e * n + j];
            const MatrixXd &B = shapes.cached_matrix_B[e * n + j];
            double det = shapes.cached_det_J[e * n + j];

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

void stiffness_stress::generate_residue(const Eigen::Ref<MatrixXd> &u, const std::vector<double> &avg_c,
    Eigen::Ref<VectorXd> res, const Eigen::Ref<Eigen::SparseMatrix<double>>& K) {
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
        if (mesh.anode_element_set.contains(e) || mesh.cathode_element_set.contains(e)) {
            VectorXd avg = VectorXd::Zero(n);
            VectorXd e_load = VectorXd::Zero(n * dim);
            auto [E, nu] = this->get_material_property(e);
            const double E_ref = 70e9;

            double omega = constant::anode.omega;
            double c_int = constant::anode.c_int;
            if (mesh.cathode_element_set.contains(e)) {
                omega = constant::cathode.omega;
                c_int = constant::cathode.c_int;
            }

            MatrixXd C(dim_voigt, dim_voigt);
            C << 1 - nu, nu, 0,
                nu, 1 - nu, 0,
                0, 0, (1 - 2 * nu) / 2;
            C = C * E / (1 + nu) / (1 - 2 * nu);
            double c_coeff = E / (1 - 2 * nu);

            for (int j = 0; j < n; j++) {
                int node_idx = dof.particle_mapper[mesh.elements[e * n + j]];
                avg(j) = avg_c[node_idx];
            }

            for (int j = 0; j < n; j++) {
                const MatrixXd &N = shapes.cached_matrix_N[e * n + j];
                const MatrixXd &B = shapes.cached_matrix_B[e * n + j];
                double det = shapes.cached_det_J[e * n + j];

                double disp = omega / 3 * ((N.transpose() * avg).value() - c_int);
                VectorXd F = VectorXd::Zero(dim_voigt);
                for (int l = 0; l < dim; l++) {
                    F(l) = disp * c_coeff;
                }

                e_load += B.transpose() * F * w(j) * det / E_ref;
            }
            //std::cout<<e_load.transpose() * E_ref<<"\n";

            for (int j = 0; j < n * dim; j++) {
                std::size_t id_l = mesh.elements[e * n + j / dim];
                std::size_t dof_l = 5 + j % dim;
                if (!mesh.anode_cc_wall_nodes.contains(id_l)) {
                    load(id_l * dim + dof_l - 5) += e_load(j);
                }
            }
        }
    }

    if(settings::calc_temperature) {
        for (int e = 0; e < mesh.elem_count; e++) {
            VectorXd e_load = VectorXd::Zero(n * dim);
            auto [E, nu] = this->get_material_property(e);
            double alpha = 1e-5;
            double c_coeff = E / (1 - 2 * nu);

            VectorXd e_t(n);
            for (int j = 0; j < n; j++) {
                std::size_t node_id = mesh.elements[e * n + j];
                e_t(j) = u(dof.get_dof(node_id, 4), 0);
            }

            for (int j = 0; j < n; j++) {
                const MatrixXd &N = shapes.cached_matrix_N[e * n + j];
                const MatrixXd &B = shapes.cached_matrix_B[e * n + j];
                double det = shapes.cached_det_J[e * n + j];
                double t = (N.transpose() * e_t).value();

                double disp = alpha * (t - constant::t_ref);
                VectorXd F = VectorXd::Zero(dim_voigt);
                for (int l = 0; l < dim; l++) {
                    F(l) = disp * c_coeff;
                }

                e_load += B.transpose() * F * w(j) * det / E_ref;
            }

            for (int j = 0; j < n * dim; j++) {
                std::size_t id_l = mesh.elements[e * n + j / dim];
                std::size_t dof_l = 5 + j % dim;
                if (!mesh.anode_cc_wall_nodes.contains(id_l)) {
                    load(id_l * dim + dof_l - 5) += e_load(j);
                }
            }
        }
    }

    /*
    for (int e = 0; e < mesh.elem_count; e++) {
        if (mesh.anode_element_set.contains(e)) {
            VectorXd avg = VectorXd::Zero(n);
            VectorXd e_load = VectorXd::Zero(n * dim);
            double E = 70e9, nu = .26, omega = -7.28e-7;
            const double E_ref = 70e9;
            double c_int = constant::c_int_an;
            if (mesh.cathode_element_set.contains(e)) {
                E = 12e9;
                nu = .3;
                omega = 4e-6;
                c_int = constant::c_int_ca;
            }

            MatrixXd C(dim_voigt, dim_voigt);
            C << 1 - nu, nu, 0,
             nu, 1 - nu, 0,
             0, 0, (1 - 2 * nu) / 2;
            C = C * E / (1 + nu) / (1 - 2 * nu);
            double c_coeff = E / (1 - 2 * nu);

            for (int j = 0; j < n; j++) {
                int node_idx = dof.particle_mapper[mesh.elements[e * n + j]];
                avg(j) = avg_c[node_idx];
            }

            for (int j = 0; j < n; j++) {
                const MatrixXd &N = shapes.cached_matrix_N[e * n + j];
                const MatrixXd &B = shapes.cached_matrix_B[e * n + j];
                double det = shapes.cached_det_J[e * n + j];

                double disp = 0.01;
                VectorXd F = VectorXd::Zero(dim_voigt);
                for (int l = 0; l < dim; l++) {
                    F(l) = disp * c_coeff;
                }
                //std::cout<<F<<"\n";
                //if (j == 0 && mesh.cathode_element_set.contains(e)) std::cout<<disp<<std::endl;

                e_load += B.transpose() * F * w(j) * det / E_ref;
            }
            //std::cout<<e_load.transpose() * E_ref<<"\n";

            for (int j = 0; j < n * dim; j++) {
                std::size_t id_l = mesh.elements[e * n + j / dim];
                std::size_t dof_l = 5 + j % dim;
                if (!mesh.anode_cc_wall_nodes.contains(id_l)) {
                    load(id_l * dim + dof_l - 5) += e_load(j);
                }
            }
        }
    }
    */

    /*
    for (int e = 0; e < mesh.elem_count; e++) {
        if (true) {
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
    }
    */

    VectorXd local_res = K * local_u - load;
    for (int x = 0; x < node_cnt; x++) {
        for (int j = 0; j < dim; j++) {
            res(dof.get_dof(x, j + 5)) = local_res(x * dim + j);
        }
    }
}

void stiffness_stress::stress_output(const Eigen::Ref<MatrixXd> &u, std::vector<double> &s) const {
    const int dim = get_dim(mesh.p_type), n = get_nodes(mesh.p_type), dim_voigt = (dim * (dim + 1)) / 2;
    std::vector<int> cnt(mesh.node_count, 0);
    s = std::vector<double>(4 * mesh.node_count, 0);

    MatrixXd T;
    if (dim == 2) {
        if (n == 3) T = get_extrapolation_weight<2, 3>();
        else if (n == 4) T = get_extrapolation_weight<2, 4>();
    }

    for (int e = 0; e < mesh.elem_count; e++) {
        auto [E, nu] = this->get_material_property(e);
        MatrixXd C(dim_voigt, dim_voigt);
        C << 1 - nu, nu, 0,
            nu, 1 - nu, 0,
            0, 0, (1 - 2 * nu) / 2;
        C = C * E / (1 + nu) / (1 - 2 * nu);

        VectorXd e_d(dim * n);

        for (int j = 0; j < n; j++) {
            int idx = mesh.elements[e * n + j];
            for (int l = 0; l < dim; l++) {
                e_d(j * dim + l) = u(dof.get_dof(idx, l + 5), 0);
            }
        }

        MatrixXd stress(dim_voigt, n);
        VectorXd strain(dim_voigt);
        VectorXd z_stress(n);

        for (int j = 0; j < n; j++) {
            const MatrixXd &B = shapes.cached_matrix_B[e * n + j];
            //double det = shapes.cached_det_J[e * n + j];

            strain = B * e_d;
            stress.col(j) = C * strain;
            z_stress(j) = E * nu / ((1 + nu) * (1 - 2 * nu)) * (strain(0) + strain(1));
        }

        MatrixXd ext_stress = (T * stress.transpose()).transpose();
        VectorXd ext_z_stress = T * z_stress;
        //std::cout<<ext_z_stress<<"\n\n";

        for (int j = 0; j < 4; j++) {
            int idx = mesh.elements[e * n + j];
            cnt[idx]++;
            for (int l = 0; l < 4; l++) {
                if (l == 0 || l == 1) {
                    s[idx * 4 + l] += ext_stress(l, j);
                } else if (l == 2) {
                    s[idx * 4 + l] += ext_z_stress(j);
                } else {
                    s[idx * 4 + l] += ext_stress(l - 1, j);
                }
            }
        }
    }

    for (int i = 0; i < mesh.node_count; i++) {
        if (cnt[i] != 0) {
            for (int j = 0; j < 4; j++) {
                s[i * 4 + j] /= cnt[i];
            }
        }
    }
}
