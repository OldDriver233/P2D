#ifndef FEM_STIFFNESS_ANODE_H
#define FEM_STIFFNESS_ANODE_H
#include "../../functions/function_manager.h"
#include "../../step_control/StepControl.h"
#include "stiffness_base.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <utility>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_anode : public stiffness_base {
public:
    double dc_ssdj;
    FunctionManager *pfm;
    StepControl* st;

    stiffness_anode() {
    }

    stiffness_anode(VectorXd points, int an, int ca, int ancoll, int cacoll, double dc_ssdj,
                    FunctionManager *pf, StepControl* st)
        : stiffness_base(std::move(points), an, ca, ancoll, cacoll), dc_ssdj(dc_ssdj), pfm(pf), st(st) {
    }

    ~stiffness_anode() {
    }

    template<bool use_temp>
    void generate(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                  bool, std::vector<double>&);
};

#endif // FEM_STIFFNESS_ANODE_H
