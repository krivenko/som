# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# start configuration 
cmake_minimum_required(VERSION 2.8)
project(som CXX)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

# Compile C++ code
add_subdirectory(c++)

# TEST only with python ?

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python)
 #add_subdirectory(test)
 option(BUILD_DOC "Build documentation" OFF)
 if(${BUILD_DOC})
   if(NOT TRIQS_WITH_DOCUMENTATION)
     message("Error: TRIQS library has not been compiled with its documentation")
   endif()
   add_subdirectory(doc)
 endif(${BUILD_DOC})
endif()

