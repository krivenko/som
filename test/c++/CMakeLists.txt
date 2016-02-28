# load triqs helper to set up tests
find_package(TriqsTest)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../c++ ${TRIQS_INCLUDE_ALL})

# Simple tests
file(GLOB test_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(filename ${test_sources})
    get_filename_component(t ${filename} NAME_WE)
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${filename})
    target_link_libraries(${t} som_c ${TRIQS_LIBRARY_ALL} gtest)
    triqs_set_rpath_for_target(${t})
    triqs_add_cpp_test(${t})
endforeach()

# MPI tests
foreach(TEST_MPI_NUMPROC 2 3 4)
    triqs_add_cpp_test(global_counter)
endforeach()

# Copy h5 files with input data
file(COPY solution_worker.ref.h5 DESTINATION .)
