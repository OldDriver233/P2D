#ifndef FEM_STIFFNESS_SEPARATOR_H
#define FEM_STIFFNESS_SEPARATOR_H
#include "../../functions/function_manager.h"
#include "../../mesh/mesh_reader.h"
#include "../../mesh/dof_assigner.h"
#include "pre_calc.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <utility>

#include "../../step_control/StepControl.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_separator {
public:
    const mesh_reader& mesh;
    const dof_assigner& dof;
    const pre_calc_shapes& shapes;
    FunctionManager *pfm;
    StepControl* st;

    stiffness_separator(const mesh_reader& mesh, const dof_assigner& dof, const pre_calc_shapes& shapes,
                        FunctionManager *pf, StepControl* st)
        : mesh(mesh), dof(dof), shapes(shapes), pfm(pf), st(st) {
    }

    ~stiffness_separator() {
    }

    template<bool use_temp>
    void generate(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                  bool, std::vector<double>&);
};

#endif // FEM_STIFFNESS_SEPARATOR_H
