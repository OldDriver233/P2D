#ifndef FEM_FUNCTION_MANAGER_H
#define FEM_FUNCTION_MANAGER_H
#include "./parser/parser.h"
#include "../io/settings/settings.h"

class FunctionManager{
public:
    Parser uoc_anode;
    Parser uoc_cathode;
    Parser anode_entropy;
    Parser cathode_entropy;
    Parser kappa;
    Parser electrolyte_diffuse;

    FunctionManager() {
        if(settings::use_customize_uoc) {
            uoc_anode.vector_size = 1;
            uoc_anode.init(settings::uoc_anode_path);
            uoc_cathode.vector_size = 1;
            uoc_cathode.init(settings::uoc_cathode_path);
            anode_entropy.vector_size = 1;
            anode_entropy.init(settings::anode_entropy_path);
            cathode_entropy.vector_size = 1;
            cathode_entropy.init(settings::cathode_entropy_path);
        }
        if (settings::use_customize_kappa) {
            kappa.vector_size = 2;
            kappa.init(settings::kappa_path);
        }
    };
};

#endif //FEM_FUNCTION_MANAGER_H