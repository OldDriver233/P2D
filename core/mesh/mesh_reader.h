#ifndef FEM_MESH_READER_H
#define FEM_MESH_READER_H 
#include <string>
#include "mesh.h"


class mesh_reader {
public:
    std::string filename;

    mesh_reader() = default;
    explicit mesh_reader(std::string filename): filename(std::move(filename)) {}

    void read(mesh&);
};


#endif // FEM_MESH_READER_H