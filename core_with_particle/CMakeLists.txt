cmake_minimum_required(VERSION 3.15...3.29)
project(fem)



#set(PYBIND11_FINDPYTHON ON)
#find_package(pybind11 CONFIG REQUIRED)

option(DEBUG "Enable debug" OFF)
option(USE_BLAS_LAPACK "Enable blas and lapack" OFF)
option(USE_MKL "Enable MKL Library" OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(IS_MACOS OFF)

if(DEBUG)
        set(CMAKE_CXX_FLAGS "-O0 -pg -Wall -Wextra -g")
else()
        set(CMAKE_CXX_FLAGS "-O3 -march=native -DNDEBUG -flto")
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        SET(IS_MACOS ON)
endif()

if (USE_BLAS_LAPACK)
        add_definitions(-DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE)
        if(NOT IS_MACOS)
                find_package(BLAS REQUIRED)
                find_package(LAPACK REQUIRED)
                find_package(OpenMP REQUIRED)
        endif()
endif()


set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(autodiff REQUIRED)

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(/usr/include/suitesparse)

add_executable(fem main.cpp
        integration/integration_shapes.h
        shaping/primitives.h
        shaping/primitive_type.h
        constants/constant.cpp
        constants/constant.h
        full_cell/full_cell_solver.cpp
        full_cell/full_cell_solver.h
        full_cell/stiffness/stiffness_separator.h
        full_cell/stiffness/stiffness_separator.cpp
        full_cell/stiffness/stiffness_anode.h
        full_cell/stiffness/stiffness_anode.cpp
        full_cell/stiffness/stiffness_cathode.h
        full_cell/stiffness/stiffness_cathode.cpp
        full_cell/stiffness/stiffness_anode_collector.cpp
        full_cell/stiffness/stiffness_cathode_collector.cpp
        full_cell/particle/particle_solver.h
        full_cell/particle/particle_solver.cpp
        functions/parser/parser.cpp
        functions/functions.h
        io/redis_connector.h
        io/redis_connector.cpp
        io/settings/settings.cpp
        io/output/output_manager.h
        io/output/output_manager.cpp
        step_control/StepControl.cpp
        mesh/mesh_reader.cpp
        mesh/mesh_reader.h
        mesh/dof_assigner.cpp
        mesh/dof_assigner.h
        full_cell/stiffness/pre_calc.h
        full_cell/apply_boundary.cpp
)

if (USE_MKL)
        add_definitions(-DEIGEN_USE_MKL_ALL -DEIGEN_DONT_PARALLELIZE -DEIGEN_UNROLLING_LIMIT=32)
        find_package(OpenMP REQUIRED)
        find_package(MKL CONFIG REQUIRED)
        find_package(SuiteSparse_config REQUIRED)
        target_compile_options(fem PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
        target_include_directories(fem PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(fem PUBLIC $<LINK_ONLY:MKL::MKL>)
        #message(ERROR ${CMAKE_CXX_FLAGS})
endif()

find_path(HIREDIS_HEADER hiredis)
target_include_directories(fem PUBLIC ${HIREDIS_HEADER})

find_library(HIREDIS_LIB hiredis)

if(NOT IS_MACOS)
target_link_libraries(fem PRIVATE nlohmann_json::nlohmann_json ${HIREDIS_LIB} "libstdc++.so" "libumfpack.so" "libblas.so")
else()
        execute_process(COMMAND /opt/homebrew/Cellar/llvm/20.1.6/bin/llvm-config --cxxflags --ldflags --system-libs --libs core
                COMMAND tr "\n" " "
                OUTPUT_VARIABLE flags)

        set(CMAKE_CXX_FLAGS "${flags} ${CMAKE_CXX_FLAGS}")
        include_directories("/opt/homebrew/Cellar/llvm/20.1.6/include")
        include_directories("/opt/homebrew/include")
        target_link_directories(fem PRIVATE "/opt/homebrew/lib")
        target_link_options(fem PRIVATE -framework Accelerate /opt/homebrew/opt/lapack/lib/liblapacke.dylib)
        target_link_libraries(fem PRIVATE nlohmann_json::nlohmann_json ${HIREDIS_LIB} "libgmsh.dylib")
endif()

install(TARGETS fem DESTINATION .)
