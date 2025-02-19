#include "constant.h"
#include "../functions/functions.h"
#include "../io/settings/settings.h"
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
double constant::sigma_an_collector;
double constant::sigma_ca_collector;
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
int constant::particle_segment;
double constant::delta_u;
double constant::R = 8.3144598;
double constant::T = 298.15;
double constant::F = 96485.33289;
double constant::t_ref = 298.15;
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
double constant::I_app = 30.0;
double constant::exchange_energy_an;
double constant::exchange_energy_ca;
double constant::diffuse_energy_an;
double constant::diffuse_energy_ca;
double constant::capacity_an_collector;
double constant::capacity_an;
double constant::capacity_sep;
double constant::capacity_ca;
double constant::capacity_ca_collector;
double constant::density_an_collector;
double constant::density_an;
double constant::density_sep;
double constant::density_ca;
double constant::density_ca_collector;
double constant::lambda_an_collector;
double constant::lambda_an;
double constant::lambda_sep;
double constant::lambda_ca;
double constant::lambda_ca_collector;


void constant::read() {
    std::ifstream f(settings::constant_path);
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
    constant::sigma_an_collector = data["sigma_an_collector"];
    constant::sigma_ca_collector = data["sigma_ca_collector"];
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
    constant::particle_segment = data["particle_segment"];
    constant::capacity_an = data["capacity_an"];
    constant::capacity_an_collector = data["capacity_an_collector"];
    constant::capacity_ca = data["capacity_ca"];
    constant::capacity_ca_collector = data["capacity_ca_collector"];
    constant::capacity_sep = data["capacity_sep"];
    constant::density_an = data["density_an"];
    constant::density_an_collector = data["density_an_collector"];
    constant::density_ca = data["density_ca"];
    constant::density_ca_collector = data["density_ca_collector"];
    constant::density_sep = data["density_sep"];
    constant::lambda_an = data["lambda_an"];
    constant::lambda_an_collector = data["lambda_an_collector"];
    constant::lambda_ca = data["lambda_ca"];
    constant::lambda_ca_collector = data["lambda_ca_collector"];
    constant::lambda_sep = data["lambda_sep"];
    constant::I_app = data["I_app"];
    constant::exchange_energy_an = data["exchange_energy_an"];
    constant::exchange_energy_ca = data["exchange_energy_ca"];
    constant::diffuse_energy_an = data["diffuse_energy_an"];
    constant::diffuse_energy_ca = data["diffuse_energy_ca"];


    constant::k_ref = kappa(ce_int);

}