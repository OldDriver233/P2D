#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "../../mesh/mesh_reader.h"
#include "../../mesh/dof_assigner.h"

using Eigen::VectorXd;

class output_manager {
public:
    const mesh_reader& mesh;
    const dof_assigner& dof;
    std::vector<double> times;
    std::vector<VectorXd> vals;

    output_manager(const mesh_reader& mesh, const dof_assigner& dof): mesh(mesh), dof(dof) {}

    void snapshot(double timestamp, const Eigen::Ref<VectorXd>& u, int dof_category);
    void snapshot_value(double timestamp, double value);
    void export_to_csv(const std::string& filename) const;
};

#endif //OUTPUT_MANAGER_H
