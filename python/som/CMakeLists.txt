# Add version file
configure_file(version.py.in version.py)

# All Python files. Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Build Python module
add_cpp2py_module(som_core)
target_link_libraries(som_core som_c)

# Install Python modules
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/som)
get_property(CPP2PY_MODULES_LIST GLOBAL PROPERTY CPP2PY_MODULES_LIST)
install(TARGETS ${CPP2PY_MODULES_LIST} DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${python_sources} ${CMAKE_CURRENT_BINARY_DIR}/version.py
        DESTINATION ${PYTHON_LIB_DEST})
