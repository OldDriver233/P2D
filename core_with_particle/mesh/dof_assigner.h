#ifndef DOF_ASSIGNER_H
#define DOF_ASSIGNER_H

#include "mesh_reader.h"
#include <vector>

class dof_assigner {
public:
    std::vector<std::size_t> dof_container;
    std::vector<std::size_t> particle_mapper;
    std::vector<std::size_t> particle_to_node;
    std::size_t dof_cnt = 0;

    dof_assigner() = default;
    ~dof_assigner() = default;

    void gen_dof(const mesh_reader& mesh);
    std::size_t get_dof(std::size_t node_id, std::size_t variable_id) const;
};

#endif //DOF_ASSIGNER_H
