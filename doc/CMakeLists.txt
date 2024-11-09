################################################################################
#
# SOM: Stochastic Optimization Method for Analytic Continuation
#
# Copyright (C) 2016-2024 Igor Krivenko
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

# Generate the Sphinx config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in
               ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

# ---------------------------------
# Top Sphinx target
# ---------------------------------
if(NOT DEFINED SPHINXBUILD_EXECUTABLE)
  find_package(Sphinx 2.1.0 REQUIRED)
endif()

# Find Sphinx extensions
find_package(SphinxRTDTheme 0.3.0 REQUIRED)
find_package(MySTParser REQUIRED)

# Sphinx has internal caching, always run it
add_custom_target(${PROJECT_NAME}_docs_sphinx ALL)
add_custom_command(
  TARGET ${PROJECT_NAME}_docs_sphinx
  COMMAND PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH}
          ${SPHINXBUILD_EXECUTABLE} -j auto -c . -b html
          ${CMAKE_CURRENT_SOURCE_DIR} html
)

# Autodoc usage requires the python modules to be built first
get_property(CPP2PY_MODULES_LIST GLOBAL PROPERTY CPP2PY_MODULES_LIST)
if(CPP2PY_MODULES_LIST)
  add_dependencies(${PROJECT_NAME}_docs_sphinx ${CPP2PY_MODULES_LIST})
endif()

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
  COMPONENT documentation
  DESTINATION share/doc/${PROJECT_NAME}
  FILES_MATCHING
  REGEX "\\.(html|pdf|png|gif|jpg|svg|js|xsl|css|py|txt|inv|bib|ttf|woff2|eot)$"
  PATTERN "_*"
)
