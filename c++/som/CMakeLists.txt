################################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2026 Igor Krivenko
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

# Find all .cpp files in subdirectories
file(GLOB_RECURSE sources *.cpp)

# Set up compilation of those files into a shared library
add_library(${PROJECT_NAME}_c ${sources})
add_library(${PROJECT_NAME}::${PROJECT_NAME}_c ALIAS ${PROJECT_NAME}_c)

# Link against triqs and enable warnings
target_link_libraries(${PROJECT_NAME}_c PUBLIC triqs nda::blas_lapack
                      PRIVATE $<BUILD_INTERFACE:${PROJECT_NAME}_warnings>)

# Configure target and compilation
set_target_properties(${PROJECT_NAME}_c PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)
target_include_directories(${PROJECT_NAME}_c PUBLIC
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
target_include_directories(${PROJECT_NAME}_c SYSTEM INTERFACE
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
target_compile_definitions(${PROJECT_NAME}_c PUBLIC
                           SOM_GIT_HASH=${PROJECT_GIT_HASH}
                           TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
                           $<$<CONFIG:Debug>:SOM_DEBUG>
                           $<$<CONFIG:Debug>:TRIQS_DEBUG>
                           $<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
                           $<$<CONFIG:Debug>:NDA_ENFORCE_BOUNDCHECK>
                          )

# Install library and headers
install(TARGETS ${PROJECT_NAME}_c
        EXPORT ${PROJECT_NAME}-targets DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
        PATTERN "CMakeFiles" EXCLUDE
)

# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES
  "Run static analyzer checks if found (clang-tidy, cppcheck)" OFF
)
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
    set_target_properties(${PROJECT_NAME}_c PROPERTIES
                          CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. "
            "Please consider installing clang-tidy for additional checks!")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
        TARGET ${PROJECT_NAME}_c
        COMMAND ${CPPCHECK_EXECUTABLE}
          --enable=warning,style,performance,portability
          --std=c++20
          --template=gcc
          --inline-suppr
          --verbose
          --force
          --quiet
          -DBOOST_PP_SEQ_FOR_EACH_PRODUCT
          ${sources}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
  else()
    message(STATUS "cppcheck not found in $PATH. "
            "Please consider installing cppcheck for additional checks!")
  endif()

endif()

# ========= Dynamic Analyzer Checks ==========

option(ASAN "Compile library with LLVM Address Sanitizer" OFF)
option(UBSAN "Compile library with LLVM Undefined Behavior Sanitizer" OFF)

if(ASAN)
  if(NOT TARGET asan)
    find_package(sanitizer REQUIRED "asan")
  endif()
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:asan>)
  install(TARGETS asan EXPORT ${PROJECT_NAME}-targets)
endif()
if(UBSAN)
  if(NOT TARGET ubsan)
    find_package(sanitizer REQUIRED "ubsan")
  endif()
  target_link_libraries(${PROJECT_NAME}_c PUBLIC ubsan)
  install(TARGETS ubsan EXPORT ${PROJECT_NAME}-targets)
endif()
