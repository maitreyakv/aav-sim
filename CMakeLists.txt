cmake_minimum_required(VERSION 3.12)
project(mpcsim)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

set(CMAKE_CXX_STANDARD 17)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(mpcsim main.cpp src/Simulation.cpp src/Simulation.h src/System.cpp src/System.h src/Dynamics.cpp src/Dynamics.h src/Controller.cpp src/Controller.h src/Cost.cpp src/Cost.h src/QuadraticCost.cpp src/QuadraticCost.h src/CartPoleDynamics.cpp src/CartPoleDynamics.h src/DDPController.cpp src/DDPController.h)

target_link_libraries(mpcsim ${ARMADILLO_LIBRARIES})
