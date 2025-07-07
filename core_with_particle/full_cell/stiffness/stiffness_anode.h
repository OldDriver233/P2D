#ifndef FEM_STIFFNESS_ANODE_H
#define FEM_STIFFNESS_ANODE_H
#include "../../functions/function_manager.h"
#include "../../step_control/StepControl.h"
#include "../../mesh/mesh_reader.h"
#include "../../mesh/dof_assigner.h"
#include "pre_calc.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <utility>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_anode {
public:
    const mesh_reader& mesh;
    const dof_assigner& dof;
    const pre_calc_shapes& shapes;
    double dc_ssdj;
    FunctionManager *pfm;
    StepControl* st;

    stiffness_anode(const mesh_reader& mesh, const dof_assigner& dof, const pre_calc_shapes& shapes, double dc_ssdj,
                    FunctionManager *pf, StepControl* st)
        : mesh(mesh), dof(dof), shapes(shapes), dc_ssdj(dc_ssdj), pfm(pf), st(st) {
    }

    ~stiffness_anode() {
    }

    template<bool use_temp>
    void generate(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  const std::vector<double>&, const std::vector<double>&,
                  std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                  bool, std::vector<double>&);
};

#endif // FEM_STIFFNESS_ANODE_H
