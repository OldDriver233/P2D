#ifndef MESH_READER_H
#define MESH_READER_H
#include <gmsh.h>
#include <Eigen/Dense>
#include <set>
#include <map>
#include "../shaping/primitive_type.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

class mesh_reader {
public:
    MatrixXd coord;
    Primitive p_type;
    std::vector<std::size_t> elements;
    std::vector<std::size_t> anode_wall;
    std::vector<std::size_t> cathode_wall;
    std::vector<std::size_t> anode_cc_elements;
    std::vector<std::size_t> cathode_cc_elements;
    std::vector<std::size_t> anode_elements;
    std::vector<std::size_t> cathode_elements;
    std::vector<std::size_t> separator_elements;
    std::set<std::size_t> anode_cc_nodes;
    std::set<std::size_t> cathode_cc_nodes;
    std::set<std::size_t> anode_nodes;
    std::set<std::size_t> cathode_nodes;
    std::set<std::size_t> separator_nodes;
    std::set<std::size_t> anode_wall_nodes;
    std::set<std::size_t> cathode_wall_nodes;
    std::vector<std::size_t> node_to_idx;
    std::size_t node_count = 0;
    std::size_t elem_count = 0;

    mesh_reader() = default;
    ~mesh_reader() = default;
    explicit mesh_reader(const std::string& filename);
};


#endif //MESH_READER_H
