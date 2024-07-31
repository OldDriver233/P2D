#ifndef FEM_CONSTANT_H
#define FEM_CONSTANT_H


class constant {
public:
    static double tolerance;
    static double dt;
    static double d_ref;
    static double d;
    static double r;
    static int step;

    static void read();
};


#endif // FEM_CONSTANT_H