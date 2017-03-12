# where will the python end up in triqs?
set(python_destination pytriqs/applications/analytical_continuation/som)

# site_customize for build
set(package_name "pytriqs.applications.analytical_continuation")
configure_file(${CMAKE_SOURCE_DIR}/cmake/sitecustomize.py ${CMAKE_CURRENT_BINARY_DIR}/sitecustomize.py @ONLY)

# make a local pytriqs copy
triqs_prepare_local_pytriqs(${python_destination})

# Build the python modules
triqs_python_extension(core ${python_destination})
target_link_libraries(core som_c ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(core)
file(COPY parameters.rst DESTINATION .)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/include ${TRIQS_INCLUDE_ALL})

install(TARGETS core DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})
install(FILES som.py __init__.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})

# add version file
configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})
