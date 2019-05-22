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

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DDEBUG_CHECK_LIMITS -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set (NETCDF_CXX "YES")
set (NETCDF_C "YES")
set (NETCDF_NC4 "YES")
FIND_PACKAGE(NetCDF)

find_package(python)
set(BLITZ_INCLUDES "lib/weave/weave/blitz")
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
set(SOURCES
  src/filters.cpp
  src/minkowski.cpp
  src/cloud_identification.cpp
  ${FILE_IO_SRC}
  src/blitz_sort.cpp
)

if(NETCDF_HAS_c++4)
  target_link_libraries(cloud_identification_core ${NETCDF_LIBRARIES})

  add_executable(main 
    ${SOURCES}
    src/main.cpp)
else()
  message(STATUS "netcdf-c++ library wasn't found, skipping compilation of CLI `main`")
endif()

add_subdirectory(lib/pybind11)
pybind11_add_module(cloud_identification
  ${SOURCE_DIR}/bindings.cpp
  ${SOURCES}
)


# For debugging cmake variables:
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#set(CMAKE_VERBOSE_MAKEFILE ON)
