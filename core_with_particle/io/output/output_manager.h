#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <eigen3/Eigen/Dense>
#include <vector>

using Eigen::VectorXd;
using Eigen::MatrixXd;

class output_manager {
public:
    std::vector<VectorXd> data;
    VectorXd coord;
    std::vector<double> tags;

    output_manager() = default;
    ~output_manager() = default;
    output_manager(const output_manager&) = default;
    explicit output_manager(VectorXd coord): coord(std::move(coord)) {}

    void append(const VectorXd& row, double tag);
    void write_to_csv(const std::string& filename);
};

#endif //OUTPUT_MANAGER_H
