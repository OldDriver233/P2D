#ifndef ACTIVE_MATERIAL_H
#define ACTIVE_MATERIAL_H
#include <fstream>
#include <nlohmann/json.hpp>

#define INPUT(entry) (this->entry = data[#entry]);

struct active_material {
    double sigma;
    double capacity;
    double density;
    double lambda;
    double E;
    double nu;
    double epsilon_e;
    double epsilon_s;
    double D_e;
    double D_s;
    double r_p;
    double c_int;
    double c_max;
    double k;
    double exchange_energy;
    double diffuse_energy;
    double omega;

    void read(const std::string& filename) {
        std::ifstream f(filename);
        nlohmann::json data = nlohmann::json::parse(f);

        INPUT(sigma)
        INPUT(capacity)
        INPUT(density)
        INPUT(lambda)
        INPUT(E)
        INPUT(nu)
        INPUT(epsilon_e)
        INPUT(epsilon_s)
        INPUT(D_e)
        INPUT(D_s)
        INPUT(r_p)
        INPUT(c_int)
        INPUT(c_max)
        INPUT(k)
        INPUT(exchange_energy)
        INPUT(diffuse_energy)
        INPUT(omega)
    }
};

#endif //ACTIVE_MATERIAL_H
