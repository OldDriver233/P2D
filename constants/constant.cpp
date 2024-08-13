#include "constant.h"
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

double constant::tolerance = 0.0;
double constant::dt = 0.0;
double constant::r = 0.0;
int constant::step = 0;
int constant::type = 0;
double constant::delta_u = 0.0;
double constant::k = 0.0;
double constant::max_c = 0.0;
double constant::initial_c = 0.0;


void constant::read() {
    std::ifstream f("config.json");
    json data = json::parse(f);
    
    tolerance = data["tolerance"];
    dt = data["dt"];
    r = data["r"];
    step = data["step"];
    delta_u = data["delta_u"];
    k = data["k"];
    max_c = data["max_c"];
    initial_c = data["initial_c"];
    type = data["type"];
}