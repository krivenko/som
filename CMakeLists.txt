# start configuration
cmake_minimum_required(VERSION 2.8)
project(triqs_som CXX)
set(CMAKE_BUILD_TYPE Release)

# Version number of the application
set(SOM_VERSION "0.9")
message(STATUS "Configuring SOM version ${SOM_VERSION}")

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

# Get hash
triqs_get_git_hash(${CMAKE_SOURCE_DIR} "SOM")
if(${GIT_RESULT} EQUAL 0)
  message(STATUS "Hash: ${SOM_GIT_HASH}")
endif(${GIT_RESULT} EQUAL 0)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

message(STATUS "TRIQS : Adding compilation flags detected by the library (C++11/14, libc++, etc...) ")
add_definitions(${TRIQS_CXX_DEFINITIONS})

# We need at least Boost 1.58 (Math Toolkit 2.2.0) for trigamma()
find_package(Boost 1.58 REQUIRED)

# Detect MPI
find_package(MPI)

set(CACHE_SIZE 0x4000 CACHE STRING "Size of the cache to store computed LHS")
add_definitions(-DCACHE_SIZE=${CACHE_SIZE})

option(EXT_DEBUG "Enable extended debugging output" OFF)
if(EXT_DEBUG)
 add_definitions(-DEXT_DEBUG)
endif(EXT_DEBUG)

option(Tests "Enable Tests" ON)

# Compile C++ code
add_subdirectory(c++)
include_directories(c++)

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python)
 if (${Tests})
  enable_testing()
  add_subdirectory(test)
 endif()
 option(BUILD_DOC "Build documentation" OFF)
 if(${BUILD_DOC})
   if(NOT ${TRIQS_WITH_DOCUMENTATION})
     message(FATAL_ERROR "Error: TRIQS library has not been compiled with its documentation")
   endif()
   add_subdirectory(doc)
 endif(${BUILD_DOC})
endif()
