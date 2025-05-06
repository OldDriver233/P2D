#include "StepControl.h"

#include <iostream>

void StepControl::step_forward() {
    if (this->mode == StepMode::Delta) {
        this->u_hist[1] = this->u_hist[0];
        this->c_s_hist[1] = this->c_s_hist[0];
    } else {
        this->u_hist[2] = this->u_hist[0];
        this->c_s_hist[2] = this->c_s_hist[0];
    }
    if (prev_time + dt_now + dt_proposed > next_output) {
        this->dt_prev = this->dt_now;
        dt_stored = dt_proposed;
        dt_now = next_output - prev_time - dt_now;
        this->status = StepStatus::ToNextOutput;
        maximum_mult = 1.0;
    } else {
        this->dt_prev = this->dt_now;
        this->dt_now = this->dt_proposed;
        this->status = StepStatus::Ok;
        maximum_mult = 2.5;
    }
    this->prev_time += this->dt_prev;
}

void StepControl::update_solution(const Eigen::Ref<VectorXd> &u, const Eigen::Ref<VectorXd>& c_s, int slot) {
    this->u_hist[slot] = u;
    this->c_s_hist[slot] = c_s;
}


void StepControl::update_dt() {
    double norm = (u_hist[1] - u_hist[0]).norm();
    double mult = 0.9 * sqrt(tolerance / norm);
    std::cout<<mult<<std::endl;
    if (mult < 0.5) {
        mult = 0.5;
        dt_now *= 0.5;
        dt_proposed = dt_now;
        maximum_mult = 1.0;
        status = StepStatus::NeedRecalc;
    } else if (status == StepStatus::NeedRecalc) {
        status = StepStatus::Ok;
    }
    mult = std::min(mult, maximum_mult);
    if (status != StepStatus::ToNextOutput) dt_proposed = dt_now * mult;
    else {
        if (mult - 1.0 > 0) dt_proposed = dt_stored;
        else dt_proposed = dt_now * mult;
    }
}
