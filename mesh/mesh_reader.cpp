#include "mesh_reader.h"
#include "mesh.h"
#include <cstddef>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <gmsh.h>
#include <utility>
#include <vector>
#include <iostream>

void gen_part(part &p, int elem_type, const std::vector<std::size_t> &elem_nodes) {
    int node_cnt;
    switch (elem_type) {
        case 1:
            p.primitive_type = Primitive::Line2;
            node_cnt = 2;
            break;
        case 2:
            p.primitive_type = Primitive::Tri3;
            node_cnt = 3;
            break;
        case 3:
            p.primitive_type = Primitive::Quad4;
            node_cnt = 4;
            break;
        case 4:
            p.primitive_type = Primitive::Tetra4;
            node_cnt = 4;
            break;
        case 5:
            p.primitive_type = Primitive::Hexa8;
            node_cnt = 8;
            break;
        case 8:
            p.primitive_type = Primitive::Line3;
            node_cnt = 3;
            break;
        default:
            break;
    }
    
    auto e = element();
    for(int i = 0; i < elem_nodes.size(); i++) {
        e.node_id.push_back(elem_nodes[i]);
        if((i + 1) % node_cnt == 0) {
            p.elements.push_back(e);
            e.node_id.clear();
        }
    }

}

void mesh_reader::read(mesh& m) {
    gmsh::initialize();
    gmsh::open(this->filename);
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);
    std::vector<std::size_t> node_tags;
    std::vector<double> node_coords, node_params;
    gmsh::model::mesh::getNodes(node_tags, node_coords, node_params);

    // The nodes might be out of order
    m.nodes = std::vector<Vector3d>(node_tags.size() + 1);
    int it = 0;
    for(auto x: node_tags) {
        m.nodes[x] = Vector3d(node_coords[3 * it], node_coords[3 * it + 1], node_coords[3 * it + 2]);
        it++;
    }

    for(auto e: entities) {
        int dim = e.first, tag = e.second;

        std::vector<int> elem_type;
        std::vector<std::vector<std::size_t>> elem_tags, elem_nodes;
        gmsh::model::mesh::getElements(elem_type, elem_tags, elem_nodes, dim, tag);

        std::vector<int> physical_tags;
        gmsh::model::getPhysicalGroupsForEntity(dim, tag, physical_tags);
        if(physical_tags.size() == 1) {
            auto phys_tag = physical_tags[0];
            std::string name;
            gmsh::model::getPhysicalName(dim, phys_tag, name);
            if(name == "left") {
                gen_part(m.bound_l, elem_type[0], elem_nodes[0]);
            } else if(name == "right") {
                gen_part(m.bound_r, elem_type[0], elem_nodes[0]);
            } else if(name == "anode") {
                gen_part(m.anode, elem_type[0], elem_nodes[0]);
            } else if(name == "cathode") {
                gen_part(m.cathode, elem_type[0], elem_nodes[0]);
            } else if(name == "separator") {
                gen_part(m.separator, elem_type[0], elem_nodes[0]);
            }
        }

    }
    gmsh::finalize();
}