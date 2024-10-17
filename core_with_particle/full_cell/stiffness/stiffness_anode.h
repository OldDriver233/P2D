#ifndef FEM_STIFFNESS_ANODE_H
#define FEM_STIFFNESS_ANODE_H
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "stiffness_base.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_anode: public stiffness_base{
public:
    double dc_ssdj;

    stiffness_anode() {}
    stiffness_anode(VectorXd points, int an, int ca, double dc_ssdj): stiffness_base(points, an, ca), dc_ssdj(dc_ssdj) {
    }
    ~stiffness_anode() {}
    void generate(const Eigen::Ref<MatrixXd>&, const Eigen::Ref<MatrixXd>&, const Eigen::Ref<MatrixXd>&, std::vector<Eigen::Triplet<double>>&, Eigen::Ref<VectorXd>, bool) override;
};

#endif //FEM_STIFFNESS_ANODE_H