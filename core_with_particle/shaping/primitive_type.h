//
// Created by paradisus on 24-7-19.
//

#ifndef FEM_PRIMITIVE_TYPE_H
#define FEM_PRIMITIVE_TYPE_H
enum class Primitive{
    Line2,
    Tri3,
    Quad4,
};

constexpr int get_dim(const Primitive p) {
    switch (p) {
        case Primitive::Line2:
            return 1;
        case Primitive::Tri3:
        case Primitive::Quad4:
            return 2;
    }
    // Unreachable
    return -1;
}

constexpr int get_nodes(const Primitive p) {
    switch (p) {
        case Primitive::Line2:
            return 2;
        case Primitive::Tri3:
            return 3;
        case Primitive::Quad4:
            return 4;
    }
    // Unreachable
    return -1;
}

#endif //FEM_PRIMITIVE_TYPE_H