#ifndef FEM_STIFFNESS_CATHODE_H
#define FEM_STIFFNESS_CATHODE_H
#include <eigen3/Eigen/Dense>
#include "stiffness_base.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class stiffness_cathode: public stiffness_base{
public:
    double dc_ssdj;

    stiffness_cathode() {}
    stiffness_cathode(VectorXd points, int an, int ca, double dc_ssdj): stiffness_base(points, an, ca), dc_ssdj(dc_ssdj) {
    }
    ~stiffness_cathode() {}
    void generate(const Eigen::Ref<MatrixXd>&, const Eigen::Ref<MatrixXd>&, const Eigen::Ref<MatrixXd>&, std::vector<Eigen::Triplet<double>>&, Eigen::Ref<VectorXd>, bool) override;
};

#endif //FEM_STIFFNESS_CATHODE_H