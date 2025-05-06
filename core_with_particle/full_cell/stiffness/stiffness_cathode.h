#ifndef FEM_STIFFNESS_CATHODE_H
#define FEM_STIFFNESS_CATHODE_H
#include "../../functions/function_manager.h"
#include "../../step_control/StepControl.h"
#include "stiffness_base.h"
#include <eigen3/Eigen/Dense>
#include <utility>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_cathode : public stiffness_base {
public:
    double dc_ssdj;
    FunctionManager *pfm;
    StepControl* st;

    stiffness_cathode() {
    }

    stiffness_cathode(VectorXd points, int an, int ca, int ancoll, int cacoll, double dc_ssdj,
                      FunctionManager *pf, StepControl* st)
        : stiffness_base(std::move(points), an, ca, ancoll, cacoll), dc_ssdj(dc_ssdj), pfm(pf), st(st) {
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
