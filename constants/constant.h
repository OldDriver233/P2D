#ifndef FEM_CONSTANT_H
#define FEM_CONSTANT_H


class constant {
public:
    static double tolerance;
    static double dt;
    static double d_ref;
    static double d;
    static double r;

    static void read() {
        tolerance = 1e-12;
        dt = 1.0;
        d_ref = 3.9e-14;
        d = 3.9e-14;
        r = 1e-5;
    }
};


#endif // FEM_CONSTANT_H