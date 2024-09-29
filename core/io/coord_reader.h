#ifndef FEM_COORD_READER_H
#define FEM_COORD_READER_H

#include <nlohmann/json.hpp>
#include <eigen3/Eigen/Dense>
#include <string>
#include <fstream>
#include <tuple>

using Eigen::VectorXd;
using json = nlohmann::json;

inline std::tuple<VectorXd, int, int> coord_reader(const std::string& filename) {
    std::ifstream f(filename);
    json data = json::parse(f);
    auto coord_vec = data["coords"].template get<std::vector<double>>();
    VectorXd coord(coord_vec.size());
    for(int i = 0; i < coord_vec.size(); i++) {
        coord(i) = coord_vec[i];
    }
    auto an = data["an"].template get<int>();
    auto ca = data["ca"].template get<int>();

    return std::make_tuple(coord, an, ca);
}

#endif //FEM_COORD_READER_H