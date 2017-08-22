cmake_minimum_required(VERSION 2.8.3)
project(cloud_identification)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -fPIC -mtune=native -funroll-all-loops -fomit-frame-pointer -march=native -mtune=native -msse4 -ftree-vectorize -ftree-vectorizer-verbose=5 -ffast-math -funroll-loops -ftracer")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fno-strict-aliasing -fwrapv -fPIC -mtune=native -funroll-all-loops -fomit-frame-pointer -march=native -mtune=native -msse4 -ftree-vectorize -ffast-math -funroll-loops -ftracer")

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release FORCE)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -DDEBUG_CHECK_LIMITS -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -Wall")

FIND_PACKAGE(PythonInterp REQUIRED)
FIND_PACKAGE(PythonLibs REQUIRED)
# I'm not sure why, but the cmake script doesn't create a variable that can be
# used when linking against the python library
set(PYTHONLIB "python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")

set (NETCDF_CXX "YES")
set (NETCDF_C "YES")
FIND_PACKAGE(NetCDF REQUIRED)

include_directories(/Library/Python/2.7/site-packages/weave/blitz/blitz ${PYTHON_INCLUDE_DIRS} ${netCDF_INCLUDE_DIR})
set(HEADER_FILES ${NETCDF_INCLUDES})

# For debugging cmake variables:
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
    #message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#set(CMAKE_VERBOSE_MAKEFILE ON)

add_executable(main cpp_identify_seungbu.cpp)
target_link_libraries(main ${PYTHONLIB} ${netCDF_LIBRARIES})
