#ifndef FEM_CONSTANT_H
#define FEM_CONSTANT_H


class constant {
public:
    static double tolerance;
    static double dt;
    static int step;
    static double epsilon_e_an;
    static double epsilon_s_an;
    static double epsilon_e_ca;
    static double epsilon_s_ca;
    static double epsilon_e_sep;
    static double epsilon_s_sep;
    static double sigma_an;
    static double sigma_ca;
    static double sigma_sep;
    static double sigma_an_collector;
    static double sigma_ca_collector;
    static double de_an;
    static double de_ca;
    static double de_sep;
    static double ds_an;
    static double ds_ca;
    static double ds_sep;
    static double r_p;
    static double l_ref;
    static double j_ref;
    static double r;
    static int type;
    static int particle_segment;
    static double delta_u;
    static double k;
    static double R, T, F;
    static double c_max_an;
    static double c_max_ca;
    static double c_int_an;
    static double c_int_ca;
    static double ce_int;
    static double k_ref;
    static double bruggeman;
    static double trans;
    static double k_an;
    static double k_ca;
    static double I_app;
    static double t_ref;
    static double exchange_energy_an;
    static double exchange_energy_ca;
    static double diffuse_energy_an;
    static double diffuse_energy_ca;
    static double capacity_an_collector;
    static double capacity_an;
    static double capacity_sep;
    static double capacity_ca;
    static double capacity_ca_collector;
    static double density_an_collector;
    static double density_an;
    static double density_sep;
    static double density_ca;
    static double density_ca_collector;
    static double lambda_an_collector;
    static double lambda_an;
    static double lambda_sep;
    static double lambda_ca;
    static double lambda_ca_collector;

    static void read();
};


#endif // FEM_CONSTANT_H