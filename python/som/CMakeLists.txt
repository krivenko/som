################################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2023 Igor Krivenko <igor.s.krivenko@gmail.com>
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

# Configure the version
configure_file(version.py.in version.py)

# Copy Python files to the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
file(GLOB_RECURSE wrap_generators RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *_desc.py)
list(REMOVE_ITEM python_sources "${wrap_generators}")
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Install Python files to proper location
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py
        DESTINATION ${PYTHON_LIB_DEST})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}
        FILES_MATCHING PATTERN "*.py" PATTERN "*_desc.py" EXCLUDE
        PATTERN "CMakeFiles" EXCLUDE)

# Build and install any python modules
foreach(gen ${wrap_generators})
  string(REPLACE "_desc.py" "" gen ${gen})
  get_filename_component(module_name ${gen} NAME_WE)
  get_filename_component(module_dir ${gen} DIRECTORY)
  add_cpp2py_module(NAME ${module_name} DIRECTORY ${module_dir})
  add_library(${PROJECT_NAME}::${module_name} ALIAS ${module_name})
  target_link_libraries(${module_name} ${PROJECT_NAME}_c triqs_py)
  install(TARGETS ${module_name} DESTINATION ${PYTHON_LIB_DEST}/${module_dir})
endforeach()
