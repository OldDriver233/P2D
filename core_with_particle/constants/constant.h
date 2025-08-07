#ifndef FEM_CONSTANT_H
#define FEM_CONSTANT_H
#include "material/active_material.h"
#include "material/current_collector.h"
#include "material/separator.h"

class constant {
public:
    static double tolerance;
    static double dt;
    static int step;
    static double output_interval;
    static double finish_time;
    static double time_step_tolerance;
    static double l_ref;
    static double j_ref;
    static double r;
    static int type;
    static int particle_segment;
    static double delta_u;
    static double k;
    static double R, T, F;
    static double ce_int;
    static double k_ref;
    static double bruggeman;
    static double trans;
    static double I_app;
    static double t_ref;
    static active_material anode, cathode;
    static current_collector anode_cc, cathode_cc;
    static separator separator;

    static void read();
};


#endif // FEM_CONSTANT_H