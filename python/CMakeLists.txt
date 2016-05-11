# make a local pytriqs copy
triqs_prepare_local_pytriqs(${CMAKE_PROJECT_NAME})

# Build the python modules
triqs_python_extension_and_install(histogram)
target_link_libraries(histogram ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(histogram)

triqs_python_extension_and_install(core)
target_link_libraries(core som_c ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(core)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/include ${TRIQS_INCLUDE_ALL})

install(FILES som.py __init__.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${CMAKE_PROJECT_NAME})

# add version file
configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${CMAKE_PROJECT_NAME})
