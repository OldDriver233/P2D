#include "constant.h"
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

double constant::tolerance = 0.0;
double constant::dt = 0.0;
double constant::d_ref = 0.0;
double constant::d = 0.0;
double constant::r = 0.0;
int constant::step = 0;

void constant::read() {
    std::ifstream f("config.json");
    json data = json::parse(f);
    
    tolerance = data["tolerance"];
    dt = data["dt"];
    d_ref = data["d_ref"];
    d = data["d"];
    r = data["r"];
    step = data["step"];
}