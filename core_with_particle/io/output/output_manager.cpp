#include "output_manager.h"

void output_manager::append(const VectorXd& row, double tag) {
    tags.push_back(tag);
    data.push_back(row);
}

void output_manager::write_to_csv(const std::string &filename) {
    std::fstream f(filename, std::ios_base::out);

    if (!f) {
        //throw std::runtime_error("Cannot open file " + filename);
    }
    f<<std::setprecision(10);

    f<<"Coord";
    for (auto i = 0; i < coord.size(); i++) {
        f<<",";
        f<<coord(i);
    }
    f<<"\n";

    int idx = 0;
    for (const auto& x: data) {
        f<<tags[idx];
        for (auto i = 0; i < x.size(); i++) {
            f<<",";
            f<<x(i);
        }
        f<<"\n";
        idx++;
    }
}
