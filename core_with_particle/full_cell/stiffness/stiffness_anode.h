#ifndef FEM_STIFFNESS_ANODE_H
#define FEM_STIFFNESS_ANODE_H
#include "../../functions/function_manager.h"
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

    stiffness_anode() {
    }

    stiffness_anode(VectorXd points, int an, int ca, int ancoll, int cacoll, double dc_ssdj,
                    FunctionManager *pf)
        : stiffness_base(std::move(points), an, ca, ancoll, cacoll), dc_ssdj(dc_ssdj), pfm(pf) {
    }

    ~stiffness_anode() {
    }

    template<bool use_temp>
    void generate(const Eigen::Ref<MatrixXd> &, const Eigen::Ref<MatrixXd> &,
                  const Eigen::Ref<MatrixXd> &,
                  std::vector<Eigen::Triplet<double> > &, Eigen::Ref<VectorXd>,
                  bool);
};

#endif // FEM_STIFFNESS_ANODE_H
