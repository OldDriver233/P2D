#ifndef FEM_STIFFNESS_CATHODE_COLLECTOR_H
#define FEM_STIFFNESS_CATHODE_COLLECTOR_H
#include "../../constants/constant.h"
#include "../../step_control/StepControl.h"
#include "stiffness_base.h"
#include <eigen3/Eigen/Dense>
#include <utility>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_cathode_collector : public stiffness_base {
public:

    StepControl* st;
    stiffness_cathode_collector() {}
    stiffness_cathode_collector(VectorXd points, int an, int ca, int ancoll,
                                int cacoll, StepControl *st)
        : stiffness_base(std::move(points), an, ca, ancoll, cacoll), st(st) {}
    ~stiffness_cathode_collector() {}

    void generate(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  std::vector<Eigen::Triplet<double>> &, Eigen::Ref<VectorXd>,
                  bool);
};

#endif
