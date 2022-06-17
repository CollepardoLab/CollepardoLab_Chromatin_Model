 # Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.13)
include(FetchContent)

# enforce out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. You must create and use a build directory. "
    "Please remove CMakeCache.txt and CMakeFiles first.")
endif()

if(NOT LAMMPS_HEADER_DIR)
  message(FATAL_ERROR "Must set LAMMPS_HEADER_DIR (e.g. cmake -DLAMMPS_HEADER_DIR=</path/to/lammps/src> <path/to/plugin/src>")
endif()

project(chromatin VERSION 1.0 LANGUAGES C CXX Fortran)

# by default, install into $HOME/.local (not /usr/local),
# so that no root access (and sudo) is needed
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE)
endif()

# C++11 is required
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Need -restrict with Intel compilers
if((CMAKE_CXX_COMPILER_ID STREQUAL "Intel") OR (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM"))
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -march=native -fno-math-errno -ffast-math -restrict")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -march=native -fno-math-errno -ffast-math")
endif()

# bail out on windows
if(CMAKE_SYSTEM_NAME STREQUAL Windows)
  message(FATAL_ERROR "LAMMPS plugins are currently not supported on Windows")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(CheckIncludeFileCXX)
include(LAMMPSInterfaceCXX)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${LAMMPS_HEADER_DIR})

# Style headers
include(StyleHeaderUtils)
RegisterStyles(${LAMMPS_HEADER_DIR})
set(LAMMPS_STYLE_HEADERS_DIR ${CMAKE_CURRENT_BINARY_DIR}/styles)
GenerateStyleHeaders(${LAMMPS_STYLE_HEADERS_DIR})
include_directories(${LAMMPS_STYLE_HEADERS_DIR})

# Set the output folder where your program will be created
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# The following folder will be included
add_subdirectory("${PROJECT_SOURCE_DIR}/lammps_custom_src")
target_include_directories(chromatin PRIVATE $ENV{HOME}/.local/include/lammps)
target_include_directories(chromatin PRIVATE "${LAMMPS_HEADER_DIR}/OPENMP")
target_link_libraries(chromatin PRIVATE lammps)

set_target_properties(chromatin PROPERTIES
                      PREFIX ""
                      LINK_FLAGS "-rdynamic")

# MacOS seems to need this
if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set_target_properties(chromatin PROPERTIES
                        LINK_FLAGS "-Wl,-undefined,dynamic_lookup")

endif()
