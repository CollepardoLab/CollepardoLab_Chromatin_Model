#===============================================================================
# @file   FindLAMMPS.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014

#
# @brief  find package for lammps
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

find_path(LAMMPS_INCLUDE_PATH lammps.h
  PATHS ${LAMMPS_DIR} ENV C_INCLUDE_PATH
  PATH_SUFFIXES src
  )

#if (not ${LAMMPS_ARCH})
  file(GLOB ARCHS "${LAMMPS_INCLUDE_PATH}/liblmp*")
  foreach(loop_var IN ITEMS ${ARCHS})
    get_filename_component(loop_var ${loop_var} NAME)
    string(REGEX REPLACE ".so" "" loop_var ${loop_var})
    string(REGEX REPLACE "liblmp_" "" loop_var ${loop_var})
#    MESSAGE ("possible archs compiled for lammps : ${loop_var}")
    SET(LAMMPS_ARCH ${loop_var} CACHE INTERNAL "internal built version of lammps detection" FORCE)
#    MESSAGE ("libname : lmp_${LAMMPS_ARCH}")
  endforeach(loop_var)
#endif(not ${LAMMPS_ARCH})


find_library(LAMMPS_MAIN_LIBRARY NAME lmp_${LAMMPS_ARCH}
  PATHS ${LAMMPS_DIR}
  PATH_SUFFIXES src
  )

if (NOT LAMMPS_MAIN_LIBRARY)
set(LAMMPS_DIR "" CACHE PATH "Location of LAMMPS library.")
endif (NOT LAMMPS_MAIN_LIBRARY)

find_library(LAMMPS_MEAM_LIBRARIES NAME meam
  PATHS ${LAMMPS_DIR}
  PATH_SUFFIXES lib/meam
)

set(LAMMPS_LIBRARIES ${LAMMPS_MAIN_LIBRARY} ${LAMMPS_MEAM_LIBRARIES})
SEPARATE_ARGUMENTS(LAMMPS_LIBRARIES)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAMMPS DEFAULT_MSG
  LAMMPS_LIBRARIES LAMMPS_INCLUDE_PATH)
