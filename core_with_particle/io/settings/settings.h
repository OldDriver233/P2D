#ifndef FEM_SETTINGS_H
#define FEM_SETTINGS_H
#include <string>

class settings {
public:
    static std::string constant_path;
    static std::string coord_path;
    static bool use_customize_uoc;
    static bool use_customize_kappa;
    static std::string uoc_anode_path;
    static std::string uoc_cathode_path;
    static std::string kappa_path;

    static void read(const std::string& filename);
};

#endif //FEM_SETTINGS_H