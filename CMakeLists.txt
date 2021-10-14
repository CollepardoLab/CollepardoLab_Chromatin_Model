#set(CMAKE_VERBOSE_MAKEFILE ON)

 # Specify the minimum version for CMake
cmake_minimum_required(VERSION 2.8)
include(FetchContent)
# Project's name

#project(chromatin VERSION 1.0 DESCRIPTION "Source code for 'Nucleosome plasticity is a critical element of chromatin liquid–liquid phase separation and multivalent nucleosome interactions'.")
project(chromatin C CXX Fortran)

#list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#find_package(LAMMPS REQUIRED)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${LAMMPS_HEADER_DIR})

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#message("${CMAKE_MODULE_PATH}")
#find_package(fmt REQUIRED)
#include_directories(${FMT_INCLUDE_PATH})

#FetchContent_Declare(fmt
#  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
#  GIT_TAG master
#)
#FetchContent_MakeAvailable(fmt)

if (${LAMMPS_FOUND})
    message("Found LAMMPS")
    message("${LAMMPS_INCLUDE_DIR}")
    message("${LAMMPS_INCLUDES}")
    message("${LAMMPS_LIBRARY}")
    message("${LAMMPS_LIBS}")
    message("${LAMMPS_LIBRARIES}")
    message("${LAMMPS_VERSION_STRING}")
endif(${LAMMPS_FOUND})

# Style headers
include(StyleHeaderUtils)
RegisterStyles(${LAMMPS_HEADER_DIR})
set(LAMMPS_STYLE_HEADERS_DIR ${CMAKE_CURRENT_BINARY_DIR}/styles)
GenerateStyleHeaders(${LAMMPS_STYLE_HEADERS_DIR})
#GeneratePackagesHeaders(${LAMMPS_STYLE_HEADERS_DIR})
include_directories(${LAMMPS_STYLE_HEADERS_DIR})

# Set the output folder where your program will be created
#set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# The following folder will be included
add_subdirectory("${PROJECT_SOURCE_DIR}/lammps_custom_src")
#include_directories("${PROJECT_SOURCE_DIR}/lammps_custom_src")

target_include_directories(chromatin PRIVATE $ENV{HOME}/.local/include/lammps)
