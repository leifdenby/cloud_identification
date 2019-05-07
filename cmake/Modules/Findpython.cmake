# This code sets the following variables:
# PYTHON_EXEC         - path to python executable
# PYTHON_LIBRARIES    - path to the python library
# PYTHON_INCLUDE_DIRS - path to where Python.h is found
# PTYHON_SITE_MODULES - path to site-packages
# PYTHON_ARCH         - name of architecture to be used for platform-specific 
#                       binary modules
# PYTHON_VERSION      - version of python
#
# You can optionally include the version number when using this package 
# like so:
#       find_package(python 2.6)
#
# This code defines a helper function find_python_module(). It can be used 
# like so:
#       find_python_module(numpy)
# If numpy is found, the variable PY_NUMPY contains the location of the numpy 
# module. If the module is required add the keyword "REQUIRED":
#       find_python_module(numpy REQUIRED)

find_program(PYTHON_EXEC "python${Python_FIND_VERSION}" 
        DOC "Location of python executable to use")
string(REGEX REPLACE "/bin/python${Python_FIND_VERSION}$" "" PYTHON_PREFIX
        "${PYTHON_EXEC}")
execute_process(COMMAND "${PYTHON_EXEC}" "-c"
        "import sys; print('%d.%d' % (sys.version_info[0],sys.version_info[1]))"
        OUTPUT_VARIABLE PYTHON_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
find_library(PYTHON_LIBRARIES "python${PYTHON_VERSION}" PATHS "${PYTHON_PREFIX}" 
        PATH_SUFFIXES "lib" "lib/python${PYTHON_VERSION}/config" 
        DOC "Python libraries" NO_DEFAULT_PATH)
find_path(PYTHON_INCLUDE_DIRS "Python.h"
        PATHS "${PYTHON_PREFIX}/include/python${PYTHON_VERSION}"
        DOC "Python include directories" NO_DEFAULT_PATH)
execute_process(COMMAND "${PYTHON_EXEC}" "-c" 
        "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
        OUTPUT_VARIABLE PYTHON_SITE_MODULES
        OUTPUT_STRIP_TRAILING_WHITESPACE)

function(find_python_module module)
        string(TOUPPER ${module} module_upper)
        if(NOT PY_${module_upper})
                if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
                        set(PY_${module}_FIND_REQUIRED TRUE)
                elseif(ARGC GREATER 1 AND ARGV1 STREQUAL "QUIET")
                        set(PY_${module}_FIND_QUIETLY TRUE)
                endif()
                # A module's location is usually a directory, but for binary modules
                # it's a .so file.
                execute_process(COMMAND "${PYTHON_EXEC}" "-c" 
                        "import re, ${module}; print(re.compile('/__init__.py.*').sub('',${module}.__file__))"
                        RESULT_VARIABLE _${module}_status 
                        OUTPUT_VARIABLE _${module}_location
                        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
                execute_process(COMMAND "${PYTHON_EXEC}" "-c" 
                        "import re, ${module}; print(${module}.__version__)"
                        OUTPUT_VARIABLE _${module}_version
                        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
                if(NOT _${module}_status)
                        set(PY_${module_upper} ${_${module}_location} CACHE STRING "Location of Python module ${module}")

                        string(REPLACE "." ";" VERSION_LIST "${_${module}_version}")
                        list(GET VERSION_LIST 0 _${module}_version_major)
                        list(GET VERSION_LIST 1 _${module}_version_minor)
                        list(GET VERSION_LIST 2 _${module}_version_patch)

                        set(PY_${module_upper}_VERSION ${_${module}_version} CACHE STRING "Version of Python module ${module}")
                        set(PY_${module_upper}_VERSION_MAJOR ${_${module}_version_major} CACHE STRING "Major version of Python module ${module}")
                        set(PY_${module_upper}_VERSION_MINOR ${_${module}_version_minor} CACHE STRING "Minor version of Python module ${module}")
                        set(PY_${module_upper}_VERSION_PATCH ${_${module}_version_patch} CACHE STRING "Patch version of Python module ${module}")
                endif(NOT _${module}_status)
        endif(NOT PY_${module_upper})
        FIND_PACKAGE_HANDLE_STANDARD_ARGS(PY_${module} DEFAULT_MSG PY_${module_upper})
endfunction(find_python_module)
