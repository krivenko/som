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

macro(FindPythonModule Name ImportName)
  execute_process(COMMAND ${CPP2PY_PYTHON_INTERPRETER}
                  "-c" "import ${ImportName}; print(${ImportName}.__version__)"
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  RESULT_VARIABLE ${Name}_FOUND
                  OUTPUT_VARIABLE ${Name}_VERSION)

  if(${${Name}_FOUND} EQUAL 0)
    set(${Name}_FOUND TRUE)
    string(REPLACE "." ";" ${Name}_VERSION_LIST ${${Name}_VERSION})
    list(GET ${Name}_VERSION_LIST 0 ${Name}_VERSION_MAJOR)
    list(GET ${Name}_VERSION_LIST 1 ${Name}_VERSION_MINOR)
    list(GET ${Name}_VERSION_LIST 2 ${Name}_VERSION_PATCH)
    unset(${Name}_VERSION_LIST)
  else(${${Name}_FOUND} EQUAL 0)
    unset(${Name}_FOUND)
    unset(${Name}_VERSION)
  endif(${${Name}_FOUND} EQUAL 0)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(${Name}
    FOUND_VAR ${Name}_FOUND
    VERSION_VAR ${Name}_VERSION
    REQUIRED_VARS ${Name}_FOUND
    FAIL_MESSAGE "Failed to find ${Name} Python module"
  )
endmacro(FindPythonModule Name ImportName)
