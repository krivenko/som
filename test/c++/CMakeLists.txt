################################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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

# Build Google Test
add_subdirectory(gtest EXCLUDE_FROM_ALL)

# Copy reference HDF5 files to binary dir
file(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
foreach(file ${all_h5_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests
file(GLOB all_tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

macro(add_cppcheck_command test_name src)
  add_custom_command(
    TARGET ${test_name}
    COMMAND ${CPPCHECK_EXECUTABLE}
      --enable=warning,style,performance,portability
      --std=c++17
      --template=gcc
      --inline-suppr
      --verbose
      --quiet
      ${src}
  )
endmacro()

# Add a test
macro(add_som_cpp_test test)
  get_filename_component(test_name ${test} NAME_WE)
  set(test_name "test.${test_name}")
  get_filename_component(test_dir ${test} DIRECTORY)
  add_executable(${test_name} ${test})

  set_property(TARGET ${test_name} PROPERTY
               RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  target_link_libraries(${test_name} som_c gtest_main project_warnings)

  if(${ARGC} GREATER 1)
    foreach(numproc ${ARGN})
      add_test(NAME ${test_name}_np${numproc}
               COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${numproc}
                       ${MPIEXEC_PREFLAGS} ${test_name} ${MPIEXEC_POSTFLAGS}
               WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
    endforeach()
  else()
    add_test(NAME ${test_name} COMMAND ${test_name}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  endif()

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    set_target_properties(${test_name} PROPERTIES
                          CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    add_cppcheck_command(${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/${test})
  endif()
endmacro()

# Simple tests
set(SIMPLE_TESTS configuration.cpp
                 objective_function.cpp
                 kernels_fermiongf.cpp
                 kernels_bosoncorr.cpp
                 kernels_bosonautocorr.cpp
                 kernels_zerotemp.cpp
                 solution_worker.cpp)
foreach(test ${SIMPLE_TESTS})
  add_som_cpp_test(${test})
endforeach()

# MPI tests
set(MPI_TESTS global_index_map.cpp spectral_stats.cpp)
foreach(test ${MPI_TESTS})
  add_som_cpp_test(${test} 1 2 4)
endforeach()

add_subdirectory(numerics)
