cmake_minimum_required(VERSION 3.26)
project(minimizer_snap_traj_generation)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -fopenmp -fPIC")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


find_package(Eigen3 REQUIRED)
find_package(OsqpEigen REQUIRED)
set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

include_directories(include
        ${Eigen3_INCLUDE_DIRS}
)

add_executable(minimizer_snap_traj_generation
        minimizer_snap/minimizer_snap_optimization.cpp
        minimizer_snap_example.cpp
)

target_link_libraries(minimizer_snap_traj_generation
        OsqpEigen::OsqpEigen
        osqp::osqp
        glog
        Python3::Python
        Python3::Module
)
