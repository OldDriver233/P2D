#ifndef FEM_STIFFNESS_SEPARATOR_H
#define FEM_STIFFNESS_SEPARATOR_H
#include "../../functions/function_manager.h"
#include "stiffness_base.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <utility>

#include "../../step_control/StepControl.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_separator : public stiffness_base {
public:
    FunctionManager *pfm;
    StepControl* st;

    stiffness_separator() {
    }

    stiffness_separator(VectorXd points, int an, int ca, int ancoll, int cacoll,
                        FunctionManager *pf, StepControl* st)
        : stiffness_base(std::move(points), an, ca, ancoll, cacoll), pfm(pf), st(st) {
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
