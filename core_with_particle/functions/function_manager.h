#ifndef FEM_FUNCTION_MANAGER_H
#define FEM_FUNCTION_MANAGER_H
#include "./parser/parser.h"
#include "../io/settings/settings.h"

class FunctionManager{
public:
    Parser uoc_anode;
    Parser uoc_cathode;
    Parser kappa;

    FunctionManager() {
        if(settings::use_customize_uoc) {
            uoc_anode.vector_size = 1;
            uoc_anode.init(settings::uoc_anode_path);
            uoc_cathode.vector_size = 1;
            uoc_cathode.init(settings::uoc_cathode_path);
        }
    };
};

#endif //FEM_FUNCTION_MANAGER_H