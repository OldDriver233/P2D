#ifndef FEM_CONSTANT_H
#define FEM_CONSTANT_H


class constant {
public:
    static double tolerance;
    static double dt;
    static double r;
    static int step;
    static int type;
    static double delta_u;
    static double k;
    static double max_c;
    static double initial_c;

    static void read();
};


#endif // FEM_CONSTANT_H