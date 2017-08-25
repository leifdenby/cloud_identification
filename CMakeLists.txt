cmake_minimum_required(VERSION 2.8.3)
project(cloud_identification)

# the version of cmake on cloud9@leeds doesn't include instructions for
# finding netcdf, so we include it here
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(BUILD_SHARED_LIBS "YES")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -fPIC -mtune=native -funroll-all-loops -fomit-frame-pointer -march=native -mtune=native -msse4 -ftree-vectorize -ftree-vectorizer-verbose=5 -ffast-math -funroll-loops -ftracer")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fno-strict-aliasing -fwrapv -fPIC -mtune=native -funroll-all-loops -fomit-frame-pointer -march=native -mtune=native -msse4 -ftree-vectorize -ffast-math -funroll-loops -ftracer")

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release FORCE)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DDEBUG_CHECK_LIMITS -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -Wall")

set (NETCDF_CXX "YES")
set (NETCDF_C "YES")
set (NETCDF_NC4 "YES")
FIND_PACKAGE(NetCDF REQUIRED)


find_package(python)
find_python_module(weave)
find_python_module(scipy.weave)

if(NOT (PY_WEAVE OR PY_SCIPY.WEAVE))
  message(FATAL_ERROR "Couldn't find python module `weave` (which is
  required to provide `blitz`). Since scipy v0.19 weave is no longer
  included with scipy. Please install with 'pip install weave'")
elseif(PY_WEAVE)
  set(BLITZ_INCLUDES ${PY_SCIPY}/weave/blitz)
elseif(PY_SCIPY.WEAVE)
  set(BLITZ_INCLUDES ${PY_SCIPY.WEAVE}/blitz)
endif()

set(SOURCE_DIR "src")


include_directories(${SOURCE_DIR} ${NETCDF_INCLUDES} ${BLITZ_INCLUDES})
set(HEADER_FILES ${NETCDF_INCLUDES})
set(SOURCES ${SOURCE_DIR})


# put all cloud identification code in library so that we link against this
# both for cli and for python module
add_library(cloud_identification
  src/cloud_identification.cpp
  src/file_io.cpp
  src/blitz_sort.cpp
)


add_executable(main src/main.cpp)
target_link_libraries(main cloud_identification ${NETCDF_LIBRARIES})


# For debugging cmake variables:
get_cmake_property(_variableNames VARIABLES)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
set(CMAKE_VERBOSE_MAKEFILE ON)
