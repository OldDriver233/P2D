#ifndef FEM_STIFFNESS_SEPARATOR_H
#define FEM_STIFFNESS_SEPARATOR_H
#include "../../functions/function_manager.h"
#include "stiffness_base.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_separator : public stiffness_base {
public:
  FunctionManager *pfm;

  stiffness_separator() {}
  stiffness_separator(VectorXd points, int an, int ca, FunctionManager *pf)
      : stiffness_base(points, an, ca), pfm(pf) {}
  ~stiffness_separator() {}
  void generate(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                const Eigen::Ref<MatrixXd> &,
                std::vector<Eigen::Triplet<double>> &, Eigen::Ref<VectorXd>,
                bool) override;
};

#endif // FEM_STIFFNESS_SEPARATOR_H