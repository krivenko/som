include(external_dependency.cmake)

if(NOT DEFINED Build_Deps)
  set(Build_Deps "Always" CACHE STRING "Do we build dependencies from source? [Never/Always/IfNotFound]")
else()
  set(Build_Deps_Opts "Never" "Always" "IfNotFound")
  if(NOT ${Build_Deps} IN_LIST Build_Deps_Opts)
    message(FATAL_ERROR "Build_Deps option should be either 'Never', 'Always' or 'IfNotFound'")
  endif()
  set(Build_Deps ${Build_Deps} CACHE STRING "Do we build dependencies from source? [Never/Always/IfNotFound]")
  if(NOT IS_SUBPROJECT AND NOT Build_Deps STREQUAL "Always" AND (ASAN OR UBSAN))
    message(WARNING "For builds with llvm sanitizers (ASAN/UBSAN) it is recommended to use -DBuild_Deps=Always to avoid false positives.")
  endif()
endif()

# -- Boost --
# We need at least Boost 1.58 (Math Toolkit 2.2.0) for trigamma()
find_package(Boost 1.58 REQUIRED)

# -- MPI --
find_package(MPI)

# -- Cpp2Py --
if(PythonSupport OR Build_Documentation)
  external_dependency(Cpp2Py
    GIT_REPO https://github.com/TRIQS/cpp2py
    VERSION 2.0
    GIT_TAG master
    BUILD_ALWAYS
    EXCLUDE_FROM_ALL
  )
endif()

# -- GTest --
external_dependency(GTest
  GIT_REPO https://github.com/google/googletest
  GIT_TAG main
  BUILD_ALWAYS
  EXCLUDE_FROM_ALL
)
