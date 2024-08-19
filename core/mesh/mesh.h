#ifndef FEM_MESH_H
#define FEM_MESH_H
#include "../shaping/primitive_type.h"
#include <cstddef>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <vector>

using Eigen::Vector3d;

struct element {
    std::vector<size_t> node_id;
};

struct part {
    Primitive primitive_type;
    std::vector<element> elements;
};

struct mesh {
    std::vector<Vector3d> nodes;
    part bound_l, bound_r, anode, cathode, separator;
};

#endif //FEM_MESH_H