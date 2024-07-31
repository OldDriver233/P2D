cmake_minimum_required(VERSION 3.15...3.29)
project(fem)

#set(PYBIND11_FINDPYTHON ON)
#find_package(pybind11 CONFIG REQUIRED)

option(DEBUG "Enable debug" OFF)

if (DEBUG)
        set(CMAKE_CXX_FLAGS "-Og -g -Wall -Wextra")
else()
        set(CMAKE_CXX_FLAGS "-O2")
endif()

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(nlohmann_json REQUIRED)

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

add_executable(fem main.cpp
        integration/integration_shapes.h
        shaping/primitives.h
        particle/particle_solver.cpp
        particle/particle_solver.h
        shaping/primitive_type.h
        particle/stiffness_generator.cpp
        particle/stiffness_generator.h
        mesh/mesh_reader.cpp
        mesh/mesh_reader.h
        constants/constant.cpp
        constants/constant.h
)

target_link_libraries(fem PRIVATE "libgmsh.so" nlohmann_json::nlohmann_json)

install(TARGETS fem DESTINATION .)
