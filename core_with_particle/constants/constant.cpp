#include "constant.h"
#include "../functions/functions.h"
#include "../io/settings/settings.h"
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

double constant::tolerance;
double constant::dt;
int constant::step;
double constant::output_interval;
double constant::finish_time;
double constant::time_step_tolerance;
double constant::l_ref;
double constant::j_ref;
double constant::r;
double constant::k;
int constant::type;
int constant::particle_segment;
double constant::delta_u;
double constant::R = 8.3144598;
double constant::T = 298.15;
double constant::F = 96485.33289;
double constant::t_ref = 298.15;
double constant::ce_int;
double constant::k_ref;
double constant::bruggeman;
double constant::trans;
double constant::I_app = 30.0;
active_material constant::anode;
active_material constant::cathode;
separator constant::separator;
current_collector constant::anode_cc;
current_collector constant::cathode_cc;

void constant::read() {
    std::ifstream f(settings::constant_path);
    json data = json::parse(f);
    
    constant::tolerance = data["tolerance"];
    constant::dt = data["dt"];
    constant::step = data["step"];
    constant::output_interval = data["output_interval"];
    if (settings::use_adaptive_time_step) {
        constant::finish_time = data["finish_time"];
    }
    constant::time_step_tolerance = data["time_step_tolerance"];
    constant::l_ref = data["l_ref"];
    constant::j_ref = data["j_ref"];
    constant::ce_int = data["ce_int"];
    constant::bruggeman = data["bruggeman"];
    constant::trans = data["trans"];
    constant::particle_segment = data["particle_segment"];
    constant::I_app = data["I_app"];
    constant::anode.read("anode.json");
    constant::cathode.read("cathode.json");
    constant::separator.read("separator.json");
    constant::anode_cc.read("anode_cc.json");
    constant::cathode_cc.read("cathode_cc.json");

    constant::k_ref = kappa(ce_int);

}