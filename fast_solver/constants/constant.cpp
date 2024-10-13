#include "constant.h"
#include "../functions/functions.h"
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

double constant::tolerance;
double constant::dt;
int constant::step;
double constant::epsilon_e_an;
double constant::epsilon_s_an;
double constant::epsilon_e_ca;
double constant::epsilon_s_ca;
double constant::epsilon_e_sep;
double constant::epsilon_s_sep;
double constant::sigma_an;
double constant::sigma_ca;
double constant::sigma_sep;
double constant::de_an;
double constant::de_ca;
double constant::de_sep;
double constant::ds_an;
double constant::ds_ca;
double constant::ds_sep;
double constant::r_p;
double constant::l_ref;
double constant::j_ref;
double constant::r;
double constant::k;
int constant::type;
double constant::delta_u;
double constant::R = 8.314;
double constant::T = 298.15;
double constant::F = 96485.3329;
double constant::c_max_an;
double constant::c_max_ca;
double constant::c_int_an;
double constant::c_int_ca;
double constant::ce_int;
double constant::k_ref;
double constant::bruggeman;
double constant::trans;
double constant::k_an;
double constant::k_ca;


void constant::read() {
    std::ifstream f("config.json");
    json data = json::parse(f);
    
    constant::tolerance = data["tolerance"];
    constant::dt = data["dt"];
    constant::step = data["step"];
    constant::epsilon_e_an = data["epsilon_e_an"];
    constant::epsilon_s_an = data["epsilon_s_an"];
    constant::epsilon_e_ca = data["epsilon_e_ca"];
    constant::epsilon_s_ca = data["epsilon_s_ca"];
    constant::epsilon_e_sep = data["epsilon_e_sep"];
    constant::epsilon_s_sep = data["epsilon_s_sep"];
    constant::sigma_an = data["sigma_an"];
    constant::sigma_ca = data["sigma_ca"];
    constant::sigma_sep = data["sigma_sep"];
    constant::de_an = data["de_an"];
    constant::de_ca = data["de_ca"];
    constant::de_sep = data["de_sep"];
    constant::ds_an = data["ds_an"];
    constant::ds_ca = data["ds_ca"];
    constant::ds_sep = data["ds_sep"];
    constant::r_p = data["r_p"];
    constant::l_ref = data["l_ref"];
    constant::j_ref = data["j_ref"];
    constant::c_max_an = data["c_max_an"];
    constant::c_max_ca = data["c_max_ca"];
    constant::c_int_an = data["c_int_an"];
    constant::c_int_ca = data["c_int_ca"];
    constant::ce_int = data["ce_int"];
    constant::k_an = data["k_an"];
    constant::k_ca = data["k_ca"];
    constant::bruggeman = data["bruggeman"];
    constant::trans = data["trans"];

    constant::k_ref = kappa(ce_int);

}