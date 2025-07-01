#ifndef FEM_STIFFNESS_ANODE_COLLECTOR_H
#define FEM_STIFFNESS_ANODE_COLLECTOR_H
#include "../../constants/constant.h"
#include "../../step_control/StepControl.h"
#include "../../mesh/mesh_reader.h"
#include "../../mesh/dof_assigner.h"
#include "pre_calc.h"
#include <Eigen/Sparse>
#include <utility>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_anode_collector {
public:
    const mesh_reader& mesh;
    const dof_assigner& dof;
    const pre_calc_shapes& shapes;
    StepControl* st;

    stiffness_anode_collector(const mesh_reader& mesh, const dof_assigner& dof, const pre_calc_shapes& shapes, StepControl* st)
        : mesh(mesh), dof(dof), shapes(shapes), st(st) {}
    ~stiffness_anode_collector() {}

    void generate(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  std::vector<Eigen::Triplet<double>> &, Eigen::Ref<VectorXd>,
                  bool);
};

#endif
