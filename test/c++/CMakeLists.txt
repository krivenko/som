# load triqs helper to set up tests
find_package(TriqsTest)

link_libraries(som_c ${TRIQS_LIBRARY_ALL} gtest)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../c++ ${TRIQS_INCLUDE_ALL})

# Simple tests
set(SIMPLE_TESTS configuration kernels objective_function)
foreach(t ${SIMPLE_TESTS})
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    triqs_set_rpath_for_target(${t})
    triqs_add_test_C_simple(${t})
endforeach(t)
