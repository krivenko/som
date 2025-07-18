################################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2025 Igor Krivenko
#
# SOM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# SOM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# SOM. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

# CMake scripts are based on the 'app4triqs' TRIQS application skeleton.
# See https://github.com/TRIQS/app4triqs for more details.

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)
cmake_policy(VERSION 3.19)
if(POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif()

# ############
# Define Project
project(som VERSION 2.1.2 LANGUAGES CXX)

get_directory_property(IS_SUBPROJECT PARENT_DIRECTORY)

# ############
# Load TRIQS
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif()
set(REQUIRED_TRIQS_VERSION_MIN 3.2)
set(REQUIRED_TRIQS_VERSION_MAX <3.4)
find_package(TRIQS
             ${REQUIRED_TRIQS_VERSION_MIN}...${REQUIRED_TRIQS_VERSION_MAX}
             REQUIRED)

# Get the git hash & print status
triqs_get_git_hash_of_source_dir(PROJECT_GIT_HASH)
message(STATUS "${PROJECT_NAME} version : ${PROJECT_VERSION}")
message(STATUS "${PROJECT_NAME} Git hash: ${PROJECT_GIT_HASH}")

# Default Install directory to TRIQS_ROOT if not given or invalid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR
   (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS
          "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT}
      CACHE PATH "Default install path" FORCE)
  set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)
endif()
if(NOT IS_SUBPROJECT)
  message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
endif()
set(${PROJECT_NAME}_BINARY_DIR ${PROJECT_BINARY_DIR}
    CACHE STRING "Binary directory of the ${PROJECT_NAME} Project")

# Make additional Find Modules available
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/share/cmake/Modules)

# ############
# CMake Options

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
if(NOT IS_SUBPROJECT)
  message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

# Python Support
option(PythonSupport "Build with Python support" ON)
if(PythonSupport AND NOT TRIQS_WITH_PYTHON_SUPPORT)
  message(FATAL_ERROR
          "TRIQS was installed without Python support. "
          "Cannot build the Python Interface. "
          "Disable the build with -DPythonSupport=OFF")
endif()

# Documentation
option(Build_Documentation "Build documentation" OFF)
if(Build_Documentation AND NOT PythonSupport)
  message(FATAL_ERROR
          "Build_Documentation=ON requires PythonSupport to be enabled")
endif()

# Testing
option(Build_Tests "Build tests" ON)
if(Build_Tests)
  enable_testing()
endif()

# Build static libraries by default
option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)

# Initial size of LHS cache
set(Cache_Size 0x4000 CACHE STRING "Size of the cache to store computed LHS")
add_definitions(-DCACHE_SIZE=${Cache_Size})

option(Ext_Debug "Enable extended debugging output" OFF)
if(Ext_Debug)
  add_definitions(-DEXT_DEBUG)
  message(STATUS "Extended debugging has been enabled")
endif(Ext_Debug)

# ############
# Global Compilation Settings

# Export the list of compile-commands into compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Disable compiler extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# Provide additional debugging information for Debug builds
add_compile_options($<$<CONFIG:Debug>:-ggdb3>)

# Create an Interface target for compiler warnings
add_library(${PROJECT_NAME}_warnings INTERFACE)
target_compile_options(${PROJECT_NAME}_warnings
  INTERFACE
    -Wall
    -Wextra
    -Wpedantic
    -Wno-sign-compare
    $<$<CXX_COMPILER_ID:GNU>:-Wno-comma-subscript>
    $<$<CXX_COMPILER_ID:GNU>:-Wshadow=local>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-attributes>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-deprecated-comma-subscript>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-unknown-warning-option>
    $<$<CXX_COMPILER_ID:Clang>:-Wshadow>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-gcc-compat>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-c++20-extensions>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wno-deprecated-comma-subscript>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wno-unknown-warning-option>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wshadow>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wno-gcc-compat>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wno-c++20-extensions>
)

# #############
# Build Project

# Find / Build dependencies
add_subdirectory(deps)

# Build and install the library
add_subdirectory(c++/${PROJECT_NAME})

# Tests
if(Build_Tests)
  add_subdirectory(test)
endif()

# Python
if(PythonSupport)
  add_subdirectory(python/${PROJECT_NAME})
endif()

# Docs
if(Build_Documentation)
  add_subdirectory(doc)
endif()

# Additional configuration files
add_subdirectory(share)

# Packaging with automatic versioning
add_subdirectory(packaging)
