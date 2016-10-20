# Numerics tests
file(GLOB test_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(filename ${test_sources})
    get_filename_component(t ${filename} NAME_WE)
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${filename})
    target_link_libraries(${t} ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST})
    triqs_set_rpath_for_target(${t})
    triqs_add_cpp_test(${t})
endforeach()
