#ifndef FEM_STIFFNESS_CATHODE_H
#define FEM_STIFFNESS_CATHODE_H
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

class stiffness_cathode {
public:
    double dc_ssdj;
    FunctionManager *pfm;
    StepControl* st;
    const pre_calc_shapes& shapes;
    const mesh_reader& mesh;
    const dof_assigner& dof;

    stiffness_cathode(const mesh_reader& mesh, const dof_assigner& dof, const pre_calc_shapes& shapes, double dc_ssdj,
                      FunctionManager *pf, StepControl* st)
        : mesh(mesh), dof(dof), dc_ssdj(dc_ssdj), pfm(pf), shapes(shapes), st(st) {
    }

    ~stiffness_cathode() {
    }

    template<bool use_temp>
    void generate(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                  bool, std::vector<double>&);
};

#endif // FEM_STIFFNESS_CATHODE_H
