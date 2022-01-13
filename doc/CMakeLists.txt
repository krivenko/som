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

# Path to the Sphinx inventory file of the library
set(TRIQS_INVENTORY_DIR "${TRIQS_ROOT}/share/doc/triqs" CACHE STRING
    "Path or URL to a directory containing a Sphinx inventory file of the TRIQS library")
mark_as_advanced(TRIQS_INVENTORY_DIR)

# generate the conf.py
configure_file(conf.py.in conf.py)

add_custom_target(docs_sphinx ALL)
add_custom_command(
  TARGET docs_sphinx
  COMMAND PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH}
          ${TRIQS_SPHINXBUILD_EXECUTABLE}
          -c . -j8 -b html ${CMAKE_CURRENT_SOURCE_DIR} html
)

# Install
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
        COMPONENT documentation
        DESTINATION share/doc/som
        REGEX "\\.(html|pdf|png|gif|jpg|svg|js|xsl|css|py|txt|inv|bib)$"
        PATTERN "_*"
)
