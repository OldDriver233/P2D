#include "settings.h"
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

std::string settings::constant_path = "config.json";
std::string settings::coord_path = "coord.json";
bool settings::use_customize_kappa = false;
bool settings::use_customize_uoc = false;
bool settings::calc_temperature = true;
bool settings::use_customize_diffuse = false;
bool settings::use_adaptive_time_step = false;
std::string settings::uoc_anode_path;
std::string settings::uoc_cathode_path;
std::string settings::kappa_path;
std::string settings::anode_entropy_path;
std::string settings::cathode_entropy_path;
std::string settings::diffuse_path;

#define TRY_READ(FIELD) if(data.contains(#FIELD)) { settings::FIELD = data[#FIELD]; }

void settings::read(const std::string& filename) {
    std::ifstream f(filename);
    if(!f.is_open()) {
        std::cout<<"No input file; using default profile"<<std::endl;
        return;
    }

    json data = json::parse(f);
    TRY_READ(use_customize_uoc)
    TRY_READ(use_customize_kappa)
    TRY_READ(use_customize_diffuse)
    TRY_READ(uoc_anode_path)
    TRY_READ(uoc_cathode_path)
    TRY_READ(kappa_path)
    TRY_READ(anode_entropy_path)
    TRY_READ(cathode_entropy_path)
    TRY_READ(diffuse_path)
    TRY_READ(calc_temperature)
    TRY_READ(use_adaptive_time_step)
}