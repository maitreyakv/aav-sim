cmake_minimum_required(VERSION 3.12)
project(quadcopter-sim)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

set(CMAKE_CXX_STANDARD 17)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

#find_package(Boost REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})

add_executable(quadcopter-sim
        src/main.cpp
        ../mpcsim/src/Simulation.cpp
        ../mpcsim/src/Simulation.h
        ../mpcsim/src/System.cpp
        ../mpcsim/src/System.h
        ../mpcsim/src/Dynamics.h
        ../mpcsim/src/Controller.h
        ../mpcsim/src/Cost.h
        ../mpcsim/src/QuadraticCost.cpp
        ../mpcsim/src/QuadraticCost.h
        ../mpcsim/src/QuadraticObstacleCost.cpp
        ../mpcsim/src/QuadraticObstacleCost.h
        ../mpcsim/src/DDPController.cpp
        ../mpcsim/src/DDPController.h
        ../mpcsim/src/SamplerController.cpp
        ../mpcsim/src/SamplerController.h
        src/QuadcopterDynamics.cpp
        src/QuadcopterDynamics.h)

target_link_libraries(quadcopter-sim ${ARMADILLO_LIBRARIES})

target_compile_options(quadcopter-sim PRIVATE -O3 -Wall -Wextra -larmadillo -DDEBUG)

if (${APPLE})
    target_compile_options(quadcopter-sim PRIVATE -framework Accelerate)
else()
    target_compile_options(quadcopter-sim PRIVATE -lblas -llapack)
endif ()
