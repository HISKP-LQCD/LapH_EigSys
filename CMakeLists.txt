# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
# Licensed under the MIT/Expat license.

# Enforce a version that is actually available on QBiG.
cmake_minimum_required(VERSION 2.8)

project(LapH-Eigsys CXX C)

include_directories("include")

# Boost
if(NOT DEFINED Boost_INCLUDE_DIRS OR NOT DEFINED Boost_LIBRARIES)
  find_package(Boost REQUIRED COMPONENTS filesystem system program_options)
endif()
include_directories(${Boost_INCLUDE_DIRS})

if(NOT DEFINED EIGEN3_INCLUDE_DIRS)
    find_package(Eigen3 REQUIRED)
endif()
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

if(NOT DEFINED PETSC_INCLUDE_DIRS OR NOT DEFINED PETSC_LIBRARIES)
   find_package(petsc REQUIRED)
endif()
include_directories(SYSTEM ${PETSC_INCLUDE_DIRS})

if(NOT DEFINED SLEPC_INCLUDE_DIRS OR NOT DEFINED SLEPC_LIBRARIES)
    find_package(slepc REQUIRED)
endif()
include_directories(SYSTEM ${SLEPC_INCLUDE_DIRS})

if(NOT DEFINED LIME_INCLUDE_DIRS OR NOT DEFINED LIME_LIBRARIES)
    find_package(lime REQUIRED)
endif()
include_directories(SYSTEM ${LIME_INCLUDE_DIRS})

find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_definitions(-Wall -pedantic --std=c++0x)
add_definitions(-Wno-unused-variable -Wno-unused-local-typedefs -Wno-sign-compare)

message(STATUS "Build Type: '${CMAKE_BUILD_TYPE}'")
if(CMAKE_BUILD_TYPE STREQUAL "")
    message(STATUS "No CMAKE_BUILD_TYPE, assuming release and enabling all optimizations.")
    set(CMAKE_BUILD_TYPE Release)
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
if(cmake_build_type_lower STREQUAL "release")
    add_definitions(-O3 -march=native)
endif()

add_executable(ev_ts
    src/ev_ts.cpp
    modules/config_utils.cpp
    modules/eigensystem.cpp
    modules/io.cpp
    modules/io_utils.cpp
    modules/ranlxs.cpp
    modules/navigation.cpp
    modules/par_io.cpp
    modules/read_write.cpp
    modules/recover_spec.cpp
    modules/shell_matop.cpp
    modules/timeslice.cpp

    )
add_executable(sourceshape
    src/sourceshape.cpp
    modules/read_write.cpp
    modules/navigation.cpp
    modules/sourceshape_funcs.cpp
    modules/par_io.cpp
    modules/timeslice.cpp
    )

add_executable(set_phase
    src/set_phase.cpp
    modules/eigensystem.cpp
    modules/navigation.cpp
    modules/read_write.cpp
    modules/sourceshape_funcs.cpp
    modules/par_io.cpp
    modules/timeslice.cpp
    )

add_executable(lapheig
    src/lapheig.cpp
    modules/io.cpp
    modules/io_utils.cpp
    modules/navigation.cpp
    modules/sourceshape_funcs.cpp
    modules/par_io.cpp
    modules/timeslice.cpp
    )

target_link_libraries(ev_ts
    ${Boost_LIBRARIES}
    ${LIME_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${SLEPC_LIBRARIES}
    m 
    )

target_link_libraries(lapheig
    ${Boost_LIBRARIES}
    ${LIME_LIBRARIES}
    m 
    )
  install(TARGETS ev_ts DESTINATION bin)
