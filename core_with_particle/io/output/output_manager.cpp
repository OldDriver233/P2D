#include "output_manager.h"
#include <fstream>
#include <iomanip>

void output_manager::snapshot(double timestamp, const Eigen::Ref<VectorXd> &u, int dof_category) {
    this->times.push_back(timestamp);
    VectorXd val(this->mesh.node_count);
    for (int i = 0; i < this->mesh.node_count; i++) {
        int dof_idx = this->dof.get_dof(i, dof_category);
        if (dof_idx == -1) {
            val(i) = NAN;
        } else {
            val(i) = u(dof_idx);
        }
    }
    this->vals.push_back(val);
}

void output_manager::snapshot_value(double timestamp, double value) {
    this->times.push_back(timestamp);
    VectorXd val(1);
    val(0) = value;
    this->vals.push_back(val);
}

void output_manager::snapshot_vector(double timestamp, const std::vector<double> &v, int dof) {
    this->times.push_back(timestamp);
    VectorXd val(this->mesh.node_count);
    for (int i = 0; i < this->mesh.node_count; i++) {
        val(i) = v[4 * i + dof];
    }
    this->vals.push_back(val);
}


void output_manager::export_to_csv(const std::string &filename) const {
    std::fstream f(filename, std::ios::out);
    f<<std::setprecision(10);
    for (int i = 0; i < this->times.size(); i++) {
        f<<times[i];
        for (int j = 0; j < this->vals[0].size(); j++) {
            f<<",";
            f<<this->vals[i][j];
        }
        f<<"\n";
    }
}
