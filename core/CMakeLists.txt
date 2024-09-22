cmake_minimum_required(VERSION 3.15...3.29)
project(fem CXX)



#set(PYBIND11_FINDPYTHON ON)
#find_package(pybind11 CONFIG REQUIRED)

option(DEBUG "Enable debug" OFF)
option(USE_BLAS_LAPACK "Enable blas and lapack" OFF)
option(USE_MKL "Enable MKL Library" OFF)

if(DEBUG)
        set(CMAKE_CXX_FLAGS "-O0 -pg -Wall -Wextra -g")
else()
        set(CMAKE_CXX_FLAGS "-O3 -march=native")
endif()

if (USE_BLAS_LAPACK)
        add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE)
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        find_package(OpenMP REQUIRED)
endif()


set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
find_package(nlohmann_json REQUIRED)

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

add_executable(fem main.cpp
        integration/integration_shapes.h
        shaping/primitives.cpp
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
        full_cell/full_cell_solver.cpp
        full_cell/full_cell_solver.h
        full_cell/stiffness/stiffness_base.h
        full_cell/stiffness/stiffness_separator.h
        full_cell/stiffness/stiffness_separator.cpp
        full_cell/stiffness/stiffness_anode.h
        full_cell/stiffness/stiffness_anode.cpp
        full_cell/stiffness/stiffness_cathode.h
        full_cell/stiffness/stiffness_cathode.cpp
        functions/functions.h
        functions/bv.cpp
        functions/j0.cpp
        functions/uoc.cpp
        functions/kappa.cpp
        io/redis_connector.h
        io/redis_connector.cpp
)

if (USE_MKL)
        add_definitions(-DEIGEN_USE_MKL_ALL)
        find_package(MKL CONFIG REQUIRED)
        target_compile_options(fem PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
        target_include_directories(fem PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(fem PUBLIC $<LINK_ONLY:MKL::MKL>)
        #message(ERROR ${CMAKE_CXX_FLAGS})
endif()

find_path(HIREDIS_HEADER hiredis)
target_include_directories(fem PUBLIC ${HIREDIS_HEADER})

find_library(HIREDIS_LIB hiredis)

target_link_libraries(fem PRIVATE "libgmsh.so" nlohmann_json::nlohmann_json ${HIREDIS_LIB} "libstdc++.so")

install(TARGETS fem DESTINATION .)
