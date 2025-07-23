#ifndef STEPCONTROL_H
#define STEPCONTROL_H
#include <eigen3/Eigen/Dense>
#include "../constants/constant.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

enum class StepStatus {
    Ok,
    NeedRecalc,
    ToNextOutput,
};

enum class StepMode {
    Delta,
    EmbeddedRK,
};

class StepControl {
public:
    double dt_stored = 0.0;
    double dt_proposed = 0.0;
    double dt_now = 1e-1;
    double dt_prev = 0.0;
    double prev_time = 0.0;
    double next_output = 0.0;
    double tolerance = constant::time_step_tolerance;
    double maximum_mult = 2.5;

    VectorXd u_hist[3];
    VectorXd c_s_hist[3];
    StepMode mode = StepMode::Delta;
    StepStatus status = StepStatus::Ok;

    StepControl() = default;
    StepControl(const StepControl&) = default;
    ~StepControl() = default;

    void step_forward();
    void update_solution(const Eigen::Ref<VectorXd>&, const Eigen::Ref<VectorXd>&, int);
    void update_dt();
};


#endif //STEPCONTROL_H
