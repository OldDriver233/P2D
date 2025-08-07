#ifndef CURRENT_COLLECTOR_H
#define CURRENT_COLLECTOR_H
#include <fstream>
#include <nlohmann/json.hpp>

#define INPUT(entry) (this->entry = data[#entry]);

struct current_collector {
    double sigma;
    double capacity;
    double density;
    double lambda;
    double E;
    double nu;

    void read(const std::string& filename) {
        std::ifstream f(filename);
        nlohmann::json data = nlohmann::json::parse(f);

        INPUT(sigma)
        INPUT(capacity)
        INPUT(density)
        INPUT(lambda)
        INPUT(E)
        INPUT(nu)
    }
};

#endif //CURRENT_COLLECTOR_H
