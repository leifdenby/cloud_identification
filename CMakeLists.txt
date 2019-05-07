cmake_minimum_required(VERSION 2.8.3)
project(cloud_identification)

# the version of cmake on cloud9@leeds doesn't include instructions for
# finding netcdf, so we include it here
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(BUILD_SHARED_LIBS "YES")

set(CMAKE_EXPORT_COMPILE_COMMANDS YES)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -fPIC -mtune=native -funroll-all-loops -fomit-frame-pointer -march=native -mtune=native -msse4 -ftree-vectorize -ftree-vectorizer-verbose=5 -ffast-math -funroll-loops -ftracer")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fno-strict-aliasing -fwrapv -fPIC -mtune=native -funroll-all-loops -fomit-frame-pointer -msse4 -ftree-vectorize -ffast-math -funroll-loops -ftracer")
# XXX: had to remove `-march=native` for compiling python interface with pybind11

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DDEBUG_CHECK_LIMITS -Wall -Warray-bounds  -ftree-vrp")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O3")

set (NETCDF_CXX "YES")
set (NETCDF_C "YES")
set (NETCDF_NC4 "YES")
FIND_PACKAGE(NetCDF)

find_package(python)
find_python_module(weave)
find_python_module(scipy.weave)

if(NOT (PY_WEAVE OR PY_SCIPY.WEAVE))
  message(FATAL_ERROR "Couldn't find python module `weave` (which is
  required to provide `blitz`). Since scipy v0.19 weave is no longer
  included with scipy. Please install with 'pip install weave'")
elseif(PY_WEAVE)
  set(BLITZ_INCLUDES ${PY_WEAVE}/blitz)
elseif(PY_SCIPY.WEAVE)
  set(BLITZ_INCLUDES ${PY_SCIPY.WEAVE}/blitz)
endif()

set(SOURCE_DIR "src")


include_directories(${SOURCE_DIR} ${NETCDF_INCLUDES} ${BLITZ_INCLUDES})

# if the netcdf c++ library isn't available we still may want to compile the
# python interface, in that case we compile in some dummy routines for loading
# data from netcdf files
if(NETCDF_HAS_c++4)
  set(FILE_IO_SRC "src/file_io.cpp")
else()
  set(FILE_IO_SRC "src/file_io_dummy.cpp")
endif()

# put all cloud identification code in library so that we link against this
# both for cli and for python module
add_library(cloud_identification_core
  src/filters.cpp
  src/minkowski.cpp
  src/cloud_identification.cpp
  ${FILE_IO_SRC}
  src/blitz_sort.cpp
)

if(NETCDF_HAS_c++4)
  target_link_libraries(cloud_identification_core ${NETCDF_LIBRARIES})

  add_executable(main 
    src/main.cpp)
  target_link_libraries(main cloud_identification_core)
else()
  message(STATUS "netcdf-c++ library wasn't found, skipping compilation of CLI `main`")
endif()

add_subdirectory(lib/pybind11)
pybind11_add_module(cloud_identification
  ${SOURCE_DIR}/bindings.cpp
)
target_link_libraries(cloud_identification cloud_identification_core)


# For debugging cmake variables:
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#set(CMAKE_VERBOSE_MAKEFILE ON)
