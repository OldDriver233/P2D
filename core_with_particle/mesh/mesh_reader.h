#ifndef MESH_READER_H
#define MESH_READER_H
#include <gmsh.h>
#include <Eigen/Dense>
#include <set>
#include "../shaping/primitive_type.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

class mesh_reader {
public:
    MatrixXd coord;
    Primitive p_type;
    std::vector<std::size_t> anode_cc;
    std::vector<std::size_t> cathode_cc;
    std::vector<std::size_t> anode;
    std::vector<std::size_t> cathode;
    std::vector<std::size_t> separator;
    std::vector<std::size_t> anode_wall;
    std::vector<std::size_t> cathode_wall;
    std::set<std::size_t> anode_cc_nodes;
    std::set<std::size_t> cathode_cc_nodes;
    std::set<std::size_t> anode_nodes;
    std::set<std::size_t> cathode_nodes;
    std::set<std::size_t> separator_nodes;
    std::set<std::size_t> anode_wall_nodes;
    std::set<std::size_t> cathode_wall_nodes;
    std::size_t node_count;

    mesh_reader() = default;
    ~mesh_reader() = default;
    explicit mesh_reader(const std::string& filename);
};


#endif //MESH_READER_H
