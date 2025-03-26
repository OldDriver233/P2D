#ifndef FEM_SETTINGS_H
#define FEM_SETTINGS_H
#include <string>

class settings {
public:
    static bool calc_temperature;
    static std::string constant_path;
    static std::string coord_path;
    static bool use_customize_uoc;
    static bool use_customize_kappa;
    static bool use_customize_diffuse;
    static bool use_adaptive_time_step;
    static std::string uoc_anode_path;
    static std::string uoc_cathode_path;
    static std::string kappa_path;
    static std::string anode_entropy_path;
    static std::string cathode_entropy_path;
    static std::string diffuse_path;

    static void read(const std::string& filename);
};

#endif //FEM_SETTINGS_H