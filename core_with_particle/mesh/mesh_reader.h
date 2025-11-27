#ifndef MESH_READER_H
#define MESH_READER_H
#include <gmsh.h>
#include <eigen3/Eigen/Dense>
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
    std::vector<std::size_t> anode_cc_wall;
    std::vector<std::size_t> cathode_cc_wall;
    std::vector<std::size_t> fixed_boundary;
    std::vector<std::size_t> x_fixed_boundary;
    std::vector<std::size_t> y_fixed_boundary;
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
    std::set<std::size_t> anode_cc_wall_nodes;
    std::set<std::size_t> cathode_cc_wall_nodes;
    std::set<std::size_t> fixed_boundary_nodes;
    std::set<std::size_t> x_fixed_boundary_nodes;
    std::set<std::size_t> y_fixed_boundary_nodes;
    std::set<std::size_t> anode_cc_element_set;
    std::set<std::size_t> cathode_cc_element_set;
    std::set<std::size_t> anode_element_set;
    std::set<std::size_t> cathode_element_set;
    std::set<std::size_t> separator_element_set;
    std::vector<std::size_t> node_to_idx;
    std::size_t node_count = 0;
    std::size_t elem_count = 0;

    mesh_reader() = default;
    ~mesh_reader() = default;
    explicit mesh_reader(const std::string& filename);
};


#endif //MESH_READER_H
