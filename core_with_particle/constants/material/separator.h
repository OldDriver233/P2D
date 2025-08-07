#ifndef SEPARATOR_H
#define SEPARATOR_H
#include <fstream>
#include <nlohmann/json.hpp>

#define INPUT(entry) (this->entry = data[#entry]);

struct separator {
    double capacity;
    double density;
    double lambda;
    double E;
    double nu;
    double epsilon_e;
    double D_e;

    void read(const std::string& filename) {
        std::ifstream f(filename);
        nlohmann::json data = nlohmann::json::parse(f);

        INPUT(capacity)
        INPUT(density)
        INPUT(lambda)
        INPUT(E)
        INPUT(nu)
        INPUT(epsilon_e)
        INPUT(D_e)
    }
};

#endif //SEPARATOR_H
