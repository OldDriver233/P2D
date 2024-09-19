#ifndef FEM_STIFFNESS_CATHODE_H
#define FEM_STIFFNESS_CATHODE_H
#include <eigen3/Eigen/Dense>
#include "stiffness_base.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_cathode: public stiffness_base{
public:
    stiffness_cathode() {}
    stiffness_cathode(VectorXd points, int an, int ca): stiffness_base(points, an, ca) {
    }
    ~stiffness_cathode() {}
    void generate(Eigen::Ref<MatrixXd>, Eigen::Ref<MatrixXd>, std::vector<Eigen::Triplet<double>>&, Eigen::Ref<VectorXd>) override;
};

#endif //FEM_STIFFNESS_CATHODE_H