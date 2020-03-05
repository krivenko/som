# Find all .cpp files in the current directory
file(GLOB cpp_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# Set up compilation of those files into a shared library
add_library(som_c STATIC ${cpp_sources})

# Link against TRIQS and enable warnings
target_link_libraries(som_c PUBLIC triqs PRIVATE project_warnings)

# Configure compilation
target_compile_options(som_c PUBLIC -fPIC)
target_include_directories(som_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
target_compile_definitions(som_c PUBLIC
  SOM_GIT_HASH=${PROJECT_GIT_HASH}
  TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
  $<$<CONFIG:Debug>:TRIQS_DEBUG>
  $<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
)

# Install library and headers
install(TARGETS som_c EXPORT som-targets DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp")

# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES OFF "Run static analyzer checks if found (clang-tidy, cppcheck)")
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
    set_target_properties(som_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
  TARGET som_c
  COMMAND ${CPPCHECK_EXECUTABLE}
    --enable=warning,style,performance,portability
    --std=c++17
    --template=gcc
    --verbose
    --force
    --quiet
    ${sources}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()


# ========= Dynamic Analyzer Checks ==========

option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")

if(ASAN)
  if(NOT TARGET asan)
    find_package(sanitizer REQUIRED "asan")
  endif()
  target_link_libraries(som_c PUBLIC asan)
  install(TARGETS asan EXPORT som-targets)
endif()
if(UBSAN)
  if(NOT TARGET ubsan)
    find_package(sanitizer REQUIRED "ubsan")
  endif()
  target_link_libraries(som_c PUBLIC ubsan)
  install(TARGETS ubsan EXPORT som-targets)
endif()
