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

# List of all tests
file(GLOB all_tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(test ${all_tests})
  add_som_cpp_test(${test})
endforeach()
