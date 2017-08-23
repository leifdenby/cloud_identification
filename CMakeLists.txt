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
FIND_PACKAGE(NetCDF REQUIRED)

find_package(python)
find_python_module(weave QUIET)
find_python_module(scipy.weave QUIET)

if(NOT (PY_WEAVE OR PY_SCIPY.WEAVE))
  message(FATAL_ERROR "Couldn't find python module `weave` (which is
  required to provide `blitz`). Since scipy v0.19 weave is no longer
  included with scipy. Please install with 'pip install weave'")
elseif(PY_WEAVE)
  set(BLITZ_INCLUDES ${PY_SCIPY}/weave/blitz)
elseif(PY_SCIPY.WEAVE)
  set(BLITZ_INCLUDES ${PY_SCIPY.WEAVE}/blitz)
endif()



include_directories(${netCDF_INCLUDE_DIR} ${BLITZ_INCLUDES})
set(HEADER_FILES ${NETCDF_INCLUDES})

# For debugging cmake variables:
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
    #message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(main cpp_identify_seungbu.cpp)
target_link_libraries(main ${NETCDF_LIBRARIES})
