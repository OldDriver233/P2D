#include "settings.h"
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

std::string settings::constant_path = "config.json";
std::string settings::coord_path = "coord.json";
bool settings::use_customize_kappa = false;
bool settings::use_customize_uoc = false;
std::string settings::uoc_anode_path = "";
std::string settings::uoc_cathode_path = "";
std::string settings::kappa_path = "";

void settings::read(const std::string& filename) {
    std::ifstream f(filename);
    if(!f.is_open()) {
        std::cout<<"No input file; using default profile"<<std::endl;
        return;
    }

    json data = json::parse(f);
    if(data.contains("constant_path")) {
        settings::constant_path = data["constant_path"];
    }
    if(data.contains("coord_path")) {
        settings::coord_path = data["coord_path"];
    }
    if(data.contains("use_customize_kappa")) {
        settings::use_customize_kappa = data["use_customize_kappa"];
    }
    if(data.contains("use_customize_uoc")) {
        settings::use_customize_uoc = data["use_customize_uoc"];
    }
    if(data.contains("uoc_anode_path")) {
        settings::uoc_anode_path = data["uoc_anode_path"];
    }
    if(data.contains("uoc_cathode_path")) {
        settings::uoc_cathode_path = data["uoc_cathode_path"];
    }
    if(data.contains("kappa_path")) {
        settings::kappa_path = data["kappa_path"];
    }
}