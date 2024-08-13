cmake_minimum_required(VERSION 3.15...3.29)
project(fem)

#set(PYBIND11_FINDPYTHON ON)
#find_package(pybind11 CONFIG REQUIRED)

option(DEBUG "Enable debug" OFF)
option(USE_BLAS_LAPACK "Enable blas and lapack" OFF)

if(DEBUG)
        set(CMAKE_CXX_FLAGS "-Og -pg -Wall -Wextra")
else()
        set(CMAKE_CXX_FLAGS "-O3 -march=native")
endif()

if (USE_BLAS_LAPACK)
        add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE)
endif()

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(OpenMP REQUIRED)

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

target_link_libraries(fem PRIVATE "libgmsh.so" nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX LAPACK::LAPACK BLAS::BLAS "liblapacke.so")

install(TARGETS fem DESTINATION .)
