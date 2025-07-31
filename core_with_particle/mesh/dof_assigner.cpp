#include "dof_assigner.h"

#include <iostream>

#include "../io/settings/settings.h"

void dof_assigner::gen_dof(const mesh_reader &mesh) {
    /// DoF generating function. The DoFs are arranged in the following order.
    /// [phi_e_1 c_e_1 phi_s_1 q_1 T_1 ... T_n]
    particle_mapper = std::vector<size_t>(mesh.node_count, -1);
    std::size_t particle_cnt = 0;
    const int dim = 2;

    dof_per_node = 5;
    if (settings::stress_analysis) {
        dof_per_node = 7;
    }
    for (int i = 0; i < mesh.node_count; i++) {
        if (mesh.anode_nodes.contains(i) || mesh.cathode_nodes.contains(i)) {
            for (int j = 0; j < 5; j++) {
                if (j == 1 || j == 4) {
                    if (j == 1 && !settings::is_solid_battery) dof_container.push_back(dof_cnt++);
                    else if (j == 4 && settings::calc_temperature) dof_container.push_back(dof_cnt++);
                    else dof_container.push_back(-1);
                } else {
                    dof_container.push_back(dof_cnt++);
                }
            }
            if (settings::stress_analysis) {
                for (int j = 0; j < dim; j++) {
                    dof_container.push_back(dof_cnt++);
                }
            }
            particle_to_node.push_back(i);
            particle_mapper[i] = particle_cnt++;
        } else if (mesh.separator_nodes.contains(i)) {
            for (int j = 0; j < 4; j++) {
                if (j == 0 || (j == 1 && !settings::is_solid_battery)) dof_container.push_back(dof_cnt++);
                else dof_container.push_back(-1);
            }
            if (settings::calc_temperature) {
                dof_container.push_back(dof_cnt++);
            } else {
                dof_container.push_back(-1);
            }
            if (settings::stress_analysis) {
                for (int j = 0; j < dim; j++) {
                    dof_container.push_back(dof_cnt++);
                }
            }
        } else {
            for (int j = 0; j < 4; j++) {
                dof_container.push_back(-1);
            }
            if (settings::calc_temperature) {
                dof_container.push_back(dof_cnt++);
            } else {
                dof_container.push_back(-1);
            }
            if (settings::stress_analysis) {
                for (int j = 0; j < dim; j++) {
                    dof_container.push_back(dof_cnt++);
                }
            }
        }
    }
}

std::size_t dof_assigner::get_dof(std::size_t node_id, std::size_t variable_id) const {
    return dof_container[dof_per_node * node_id + variable_id];
}
