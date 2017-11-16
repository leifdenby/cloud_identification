################################################################################
#
# \file      cmake/FindLNetCDF.cmake
# \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
#            2017, Leif Denby, University of Leeds
# \brief     Find NetCDF
#
################################################################################

# - Find NetCDF
# Find the native NetCDF includes and library
#
#  NETCDF_INCLUDES     - where to find netcdf.h, etc
#  NETCDF_LIBRARIES    - Link these libraries when using NetCDF
#  NETCDF_LIBRARY_DIRS - location of the NetCDF libraries
#  NETCDF_FOUND        - True if NetCDF found including required interfaces (see below)
#
# Your package can require certain interfaces to be FOUND by setting these
#
#  NETCDF_CXX         - require the C++ interface and link the C++ library
#  NETCDF_F77         - require the F77 interface and link the fortran library
#  NETCDF_F90         - require the F90 interface and link the fortran library
#
# The following are not for general use and are included in
# NETCDF_LIBRARIES if the corresponding option above is set.
#
#  NETCDF_LIBRARIES_C    - Just the C interface
#  NETCDF_LIBRARIES_CXX  - C++ interface, if available
#  NETCDF_LIBRARIES_F77  - Fortran 77 interface, if available
#  NETCDF_LIBRARIES_F90  - Fortran 90 interface, if available
#
# Normal usage would be:
#  set (NETCDF_ROOT ${CUSTOM_PATH})
#  set (NETCDF_F90 "YES")
#  find_package (NetCDF REQUIRED)
#  target_link_libraries (uses_f90_interface ${NETCDF_LIBRARIES})
#  target_link_libraries (only_uses_c_interface ${NETCDF_LIBRARIES_C})

# If already in cache, be silent
if (NETCDF_INCLUDES AND NETCDF_LIBRARIES)
  set (NETCDF_FIND_QUIETLY TRUE)
endif (NETCDF_INCLUDES AND NETCDF_LIBRARIES)

if(NETCDF_NC4)
  macro(NetCDF_check_feature feature_name)
    execute_process(COMMAND "nc-config" "--has-${feature_name}" OUTPUT_VARIABLE NETCDF_HAS_${feature_name} OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NC_CONFIG_MISSING)
      message(FATAL_ERROR "Couldn't find nc-config to determine netcdf capabilities")
    endif()
  endmacro (NetCDF_check_feature)
  NetCDF_check_feature("nc4")

  if(NOT NETCDF_HAS_nc4)
    message(FATAL_ERROR "NetCDF not compiled with v4 support")
  else()
    message(STATUS "NetCDF4 found")
  endif()

  macro(NetCDF_get_flags flags_name output_var)
    execute_process(COMMAND "nc-config" "--${flags_name}" OUTPUT_VARIABLE NETCDF_${output_var} OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_VARIABLE NETCDF_BOO)
    # annoyingly there's a version of fortran netcdf where `nc-config` doesn't work if compiled with cmake
    if(NETCDF_${output_var} STREQUAL "nf-config not yet implemented for cmake builds")
      set(NETCDF_${output_var} "")
    endif()
  endmacro (NetCDF_get_flags)

  if(NETCDF_CXX)
    NetCDF_check_feature("c++4")
    if(NOT NETCDF_HAS_c++4)
      message(STATUS "NetCDF not compiled with c++ support")
    else()
      NetCDF_get_flags("cflags" INCLUDES_CXX)
      NetCDF_get_flags("libs" LIBRARIES_CXX)
      if (NOT NETCDF_LIBRARIES_CXX)
        set(NETCDF_HAS_c++4 "NO")
      elseif (NOT NETCDF_INCLUDES_CXX)
        set(NETCDF_HAS_c++4 "NO")
      else()
        message(STATUS "Found NetCDF4 with C++ bindings")
        set (NETCDF_LIBRARIES "${NETCDF_LIBRARIES_CXX}" CACHE STRING "All NetCDF libraries required for interface level")
        set (NETCDF_INCLUDES "${NETCDF_INCLUDES_CXX}" CACHE STRING "All NetCDF libraries required for interface level")
      endif()
    endif()
  endif()

  include (FindPackageHandleStandardArgs)
  find_package_handle_standard_args (NetCDF DEFAULT_MSG NETCDF_LIBRARIES NETCDF_INCLUDES)

  mark_as_advanced (NETCDF_LIBRARIES NETCDF_INCLUDES)
else()

  find_path (NETCDF_INCLUDES netcdf.h
    HINTS ${NETCDF_ROOT}/include
    ${NETCDF_DIR}/include
    $ENV{NETCDF_DIR}/include
    $ENV{NETCDF_HOME}/include)

  if(NOT BUILD_SHARED_LIBS)
    find_library (NETCDF_LIBRARIES_C NAMES libnetcdf.a
      HINTS ${NETCDF_ROOT}/lib
      ${NETCDF_DIR}/lib
      $ENV{NETCDF_DIR}/lib
      $ENV{NETCDF_HOME}/lib)
  else()
    find_library (NETCDF_LIBRARIES_C NAMES netcdf
      HINTS ${NETCDF_ROOT}/lib
      ${NETCDF_DIR}/lib
      $ENV{NETCDF_DIR}/lib
      $ENV{NETCDF_HOME}/lib)
  endif()
  mark_as_advanced(NETCDF_LIBRARIES_C)

  set (NetCDF_has_interfaces "YES") # will be set to NO if we're missing any interfaces
  set (NetCDF_libs "${NETCDF_LIBRARIES_C}")

  get_filename_component (NETCDF_LIBRARY_DIRS "${NETCDF_LIBRARIES_C}" PATH)
  mark_as_advanced(NETCDF_LIBRARY_DIRS)

  macro (NetCDF_check_interface lang header libs)
    if (NETCDF_${lang})
      find_path (NETCDF_INCLUDES_${lang} NAMES ${header}
        HINTS "${NETCDF_INCLUDES}" NO_DEFAULT_PATH)
      find_library (NETCDF_LIBRARIES_${lang} NAMES ${libs}
        HINTS "${NETCDF_LIBRARY_DIRS}" NO_DEFAULT_PATH)

      mark_as_advanced (NETCDF_INCLUDES_${lang} NETCDF_LIBRARIES_${lang})
      if (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
        list (INSERT NetCDF_libs 0 ${NETCDF_LIBRARIES_${lang}}) # prepend so that -lnetcdf is last
      else (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
        set (NetCDF_has_interfaces "NO")
        message (STATUS "Failed to find NetCDF interface for ${lang}")
      endif (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
    endif (NETCDF_${lang})
  endmacro (NetCDF_check_interface)

  if(NOT BUILD_SHARED_LIBS)
    NetCDF_check_interface (CXX netcdfcpp.h libnetcdf_c++.a)
    NetCDF_check_interface (F77 netcdf.inc  libnetcdff.a)
    NetCDF_check_interface (F90 netcdf.mod  libnetcdff.a)
  else()
    NetCDF_check_interface (CXX netcdfcpp.h netcdf_c++)
    NetCDF_check_interface (F77 netcdf.inc  netcdff)
    NetCDF_check_interface (F90 netcdf.mod  netcdff)
  endif()

  set (NETCDF_LIBRARIES "${NetCDF_libs}" CACHE STRING "All NetCDF libraries required for interface level")

  # handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
  # all listed variables are TRUE
  include (FindPackageHandleStandardArgs)
  find_package_handle_standard_args (NetCDF DEFAULT_MSG NETCDF_LIBRARIES NETCDF_INCLUDES NetCDF_has_interfaces)

  mark_as_advanced (NETCDF_LIBRARIES NETCDF_INCLUDES)
endif()
