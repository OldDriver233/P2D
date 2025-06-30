#include "dof_assigner.h"
#include "../io/settings/settings.h"

void dof_assigner::gen_dof(const mesh_reader &mesh) {
    /// DoF generating function. The DoFs are arranged in the following order.
    /// [phi_e_1 c_e_1 phi_s_1 q_1 T_1 ... T_n]
    for (int i = 0; i < mesh.node_count; i++) {
        if (mesh.anode_nodes.contains(i) || mesh.cathode_nodes.contains(i)) {
            for (int j = 0; j < 4; j++) {
                dof_container.push_back(dof_cnt++);
            }
            if (settings::calc_temperature) {
                dof_container.push_back(dof_cnt++);
            }
        } else if (mesh.separator_nodes.contains(i)) {
            for (int j = 0; j < 4; j++) {
                if (j < 2) dof_container.push_back(dof_cnt++);
                else dof_container.push_back(-1);
            }
            if (settings::calc_temperature) {
                dof_container.push_back(dof_cnt++);
            }
        } else {
            for (int j = 0; j < 4; j++) {
                dof_container.push_back(-1);
            }
            if (settings::calc_temperature) {
                dof_container.push_back(dof_cnt++);
            }
        }
    }
}
