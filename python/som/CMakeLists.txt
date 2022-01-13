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

# Add version file
configure_file(version.py.in version.py)

# Copy Python files to the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
list(REMOVE_ITEM python_sources "som_core_desc.py")
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Build Python module
add_cpp2py_module(som_core)
target_link_libraries(som_core som_c triqs_py)

# Install Python modules
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/som)
get_property(CPP2PY_MODULES_LIST GLOBAL PROPERTY CPP2PY_MODULES_LIST)
install(TARGETS ${CPP2PY_MODULES_LIST} DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${python_sources} ${CMAKE_CURRENT_BINARY_DIR}/version.py
        DESTINATION ${PYTHON_LIB_DEST})
