#
# Compile the cython module and install the python files

# where will the python end up in triqs?
set(python_destination applications/my_project)

# make a local pytriqs copy to be able to run in the current dir without installation
triqs_prepare_local_pytriqs_merged_with_my_python(${python_destination})

# Build the cython module
triqs_cython_module(my_project my_project ${python_destination})

# The cython module needs to linked to the C++ lib
target_link_libraries(my_project my_project_c)

# The cython module require some include path, the same as TRIQS plus ./
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_SOURCE_DIR})

# This macro sets up the rpath, to avoid any LD_LIBRARY_PATH issue ...
triqs_set_rpath_for_target(my_project)

# to be able to run from toplevel, we add a link
execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/my_project.so ${CMAKE_BINARY_DIR}/pytriqs/${python_destination})

