#ifndef FEM_STIFFNESS_SEPARATOR_H
#define FEM_STIFFNESS_SEPARATOR_H
#include "../../functions/function_manager.h"
#include "stiffness_base.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <utility>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_separator : public stiffness_base {
public:
    FunctionManager *pfm;

    stiffness_separator() {
    }

    stiffness_separator(VectorXd points, int an, int ca, int ancoll, int cacoll,
                        FunctionManager *pf)
        : stiffness_base(std::move(points), an, ca, ancoll, cacoll), pfm(pf) {
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
