#include "mesh_reader.h"
#include <iostream>
#include <algorithm>
#include <ranges>

#define GET_ELEMENT(ty) \
gmsh::model::getEntitiesForPhysicalName(#ty, tags);\
dim = tags[0].first, tag = tags[0].second;\
gmsh::model::mesh::getElements(e_type, elem_tag, node_tag, dim, tag);\
assert(e_type.size() == 1);\
for (auto &x: node_tag[0]) {\
    x -= 1;\
    this->ty##_nodes.emplace(x);\
    this->elements.emplace_back(x);\
}\
for (auto i = 0; i < elem_tag[0].size(); i++){\
    this->ty##_elements.emplace_back(elem_count++);\
}

mesh_reader::mesh_reader(const std::string &filename) {
    gmsh::initialize();
    gmsh::open(filename);

    gmsh::model::mesh::getMaxNodeTag(node_count);

    std::vector<std::pair<int, int>> tags;
    int dim, tag;
    std::vector<int> e_type;
    std::vector<std::vector<size_t>> elem_tag;
    std::vector<std::vector<size_t>> node_tag;

    GET_ELEMENT(separator)

    switch (e_type[0]) {
        case 1: this->p_type = Primitive::Line2; this->coord = MatrixXd::Zero(1, node_count); break;
        case 2: this->p_type = Primitive::Tri3; this->coord = MatrixXd::Zero(2, node_count); break;
        case 3: this->p_type = Primitive::Quad4; this->coord = MatrixXd::Zero(2, node_count); break;
        default:
            std::cerr<<"Unimplemented shape\n";
    }

    GET_ELEMENT(anode)
    GET_ELEMENT(cathode)
    GET_ELEMENT(anode_cc)
    GET_ELEMENT(cathode_cc)
    //GET_ELEMENT(anode_wall)
    //GET_ELEMENT(cathode_wall)
    gmsh::model::getEntitiesForPhysicalName("anode_wall", tags);
    dim = tags[0].first, tag = tags[0].second;
    gmsh::model::mesh::getElements(e_type, elem_tag, node_tag, dim, tag);\
    assert(e_type.size() == 1);
    this->anode_wall = node_tag[0];
    for (auto &x: this->anode_wall) {
        x -= 1;
        this->anode_wall_nodes.emplace(x);
    }
    gmsh::model::getEntitiesForPhysicalName("cathode_wall", tags);
    dim = tags[0].first, tag = tags[0].second;
    gmsh::model::mesh::getElements(e_type, elem_tag, node_tag, dim, tag);\
    assert(e_type.size() == 1);
    this->cathode_wall = node_tag[0];
    for (auto &x: this->cathode_wall) {
        x -= 1;
        this->cathode_wall_nodes.emplace(x);
    }

    int idx = 0;
    this->node_to_idx = std::vector<std::size_t>(node_count, -1);
    for (auto x: this->anode_nodes) {
        this->node_to_idx[x] = idx++;
    }
    idx = 0;
    for (auto x: this->cathode_nodes) {
        this->node_to_idx[x] = idx++;
    }

    std::vector<size_t> node_tags;
    std::vector<double> packed_coord, local_coord;

    gmsh::model::mesh::getNodes(node_tags, packed_coord, local_coord, -1, -1, true, false);

    for (int i = 0; i < node_tags.size(); i++) {
        coord(0, node_tags[i] - 1) = packed_coord[3 * i];
        if (coord.rows() >= 2) coord(1, node_tags[i] - 1) = packed_coord[3 * i + 1];
    }

    gmsh::finalize();
}
