find_package(TriqsTest)

file(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${all_h5_files} DESTINATION .)

triqs_add_python_test(gf_imtime)
triqs_add_python_test(gf_imfreq)
triqs_add_python_test(gf_legendre)
