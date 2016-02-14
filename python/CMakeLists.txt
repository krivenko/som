# make a local pytriqs copy
triqs_prepare_local_pytriqs(${CMAKE_PROJECT_NAME})

# Build the python module
triqs_python_extension_and_install(som)
target_link_libraries(som som_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(som)

# to be able to run from toplevel
execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/som.so ${CMAKE_BINARY_DIR}/pytriqs/${python_destination})

# add version file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.py.in ${CMAKE_BINARY_DIR}/pytriqs/${python_destination}/version.py)
install (FILES ${CMAKE_BINARY_DIR}/pytriqs/${python_destination}/version.py DESTINATION ${TRIQS_PYTHON_LIB_DEST}/${python_destination})
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/${python_destination})