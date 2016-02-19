# load triqs helper to set up tests
find_package(TriqsTest)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../c++ ${TRIQS_INCLUDE_ALL})

# Simple tests
set(TESTS configuration kernels objective_function global_counter gtsv spline pool)
set(TESTS gf ${TESTS})
foreach(t ${TESTS})
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    target_link_libraries(${t} som_c ${TRIQS_LIBRARY_ALL} gtest)
    triqs_set_rpath_for_target(${t})
    triqs_add_cpp_test(${t})
endforeach(t)
