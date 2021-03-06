cmake_minimum_required(VERSION 3.12)
project(aav-sim)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

set(CMAKE_CXX_STANDARD 17)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(LEMON DEFAULT_MSG LEMON_LIBRARY LEMON_INCLUDE_DIR)

#find_package(Boost REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})

add_executable(aav-sim
        src/main.cpp
        src/mpc/control/Controller.h
        src/mpc/control/DDPController.cpp
        src/mpc/control/DDPController.h
        src/mpc/cost/Cost.h
        src/mpc/cost/QuadraticCost.cpp
        src/mpc/cost/QuadraticCost.h
        src/quadcopter/QuadcopterDynamics.cpp
        src/quadcopter/QuadcopterDynamics.h
        src/simulation/Dynamics.h
        src/simulation/Simulation.cpp
        src/simulation/Simulation.h
        src/simulation/System.cpp
        src/simulation/System.h
        src/navigation/PathPlanner.h
        src/simulation/Obstacle.h
        src/navigation/RRTPathPlanner.cpp
        src/navigation/RRTPathPlanner.h
        src/simulation/SphereObstacle.cpp
        src/simulation/SphereObstacle.h)

target_link_libraries(aav-sim ${ARMADILLO_LIBRARIES})

target_compile_options(aav-sim PRIVATE -O3 -Wall -Wextra -larmadillo -lemon -DDEBUG)

if (${APPLE})
    target_compile_options(aav-sim PRIVATE -framework Accelerate)
else()
    target_compile_options(aav-sim PRIVATE -lblas -llapack)
endif ()
