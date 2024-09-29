#ifndef FEM_STIFFNESS_SEPARATOR_H
#define FEM_STIFFNESS_SEPARATOR_H
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "stiffness_base.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_separator: public stiffness_base{
public:
    stiffness_separator() {}
    stiffness_separator(VectorXd points, int an, int ca): stiffness_base(points, an, ca) {
    }
    ~stiffness_separator() {}
    void generate(Eigen::Ref<MatrixXd>, Eigen::Ref<MatrixXd>, std::vector<Eigen::Triplet<double>>&, Eigen::Ref<VectorXd>, bool) override;
};

#endif //FEM_STIFFNESS_SEPARATOR_H