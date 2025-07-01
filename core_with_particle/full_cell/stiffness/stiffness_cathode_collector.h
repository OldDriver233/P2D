#ifndef FEM_STIFFNESS_CATHODE_COLLECTOR_H
#define FEM_STIFFNESS_CATHODE_COLLECTOR_H
#include "../../constants/constant.h"
#include "../../step_control/StepControl.h"
#include "../../mesh/mesh_reader.h"
#include "../../mesh/dof_assigner.h"
#include "pre_calc.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <utility>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_cathode_collector {
public:
    const pre_calc_shapes& shapes;
    const mesh_reader& mesh;
    const dof_assigner& dof;
    StepControl* st;
    stiffness_cathode_collector(const mesh_reader& mesh, const dof_assigner& dof, const pre_calc_shapes& shapes, StepControl *st)
        : mesh(mesh), dof(dof), shapes(shapes), st(st) {}
    ~stiffness_cathode_collector() {}

    void generate(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  std::vector<Eigen::Triplet<double>> &, Eigen::Ref<VectorXd>,
                  bool);
};

#endif
