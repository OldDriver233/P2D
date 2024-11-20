#ifndef FEM_STIFFNESS_ANODE_COLLECTOR_H
#define FEM_STIFFNESS_ANODE_COLLECTOR_H
#include "../../constants/constant.h"
#include "stiffness_base.h"
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_anode_collector : public stiffness_base {
public:
  stiffness_anode_collector() {}
  stiffness_anode_collector(VectorXd points, int an, int ca, int ancoll,
                            int cacoll)
      : stiffness_base(points, an, ca, ancoll, cacoll) {}
  ~stiffness_anode_collector() {}

  void generate(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                const Eigen::Ref<MatrixXd> &,
                std::vector<Eigen::Triplet<double>> &, Eigen::Ref<VectorXd>,
                bool) override;
};

#endif
