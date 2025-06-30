#ifndef DOF_ASSIGNER_H
#define DOF_ASSIGNER_H

#include "mesh_reader.h"
#include <vector>

class dof_assigner {
public:
    std::vector<std::size_t> dof_container;
    std::size_t dof_cnt = 0;

    dof_assigner() = default;
    ~dof_assigner() = default;

    void gen_dof(const mesh_reader& mesh);
};

#endif //DOF_ASSIGNER_H
