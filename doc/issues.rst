.. _issues:

Reporting issues
================

Please report all problems and bugs directly at the GitHub issue page
`<https://github.com/krivenko/som/issues>`_.  In order to make it easier for me
to solve the issue please follow these guidelines:

#. In all cases specify which version of the application you are using. You can
   find the version number in the file ``CMakeLists.txt`` at the root of the
   application sources. It is also reported by CMake early in the configuration
   process. If the application is already installed, you can learn its version
   number by running

   ::

      $ pytriqs -c "from pytriqs.applications.analytical_continuation.som import version; version.show_version()"

#. If you have a problem during the installation, give me information about
   your operating system and the compiler you are using. Include the outputs of
   the ``cmake`` and ``make`` commands as well as the ``CMakeCache.txt`` file
   which is in the build directory. Please include these outputs in a
   `gist <http://gist.github.com/>`_ file referenced in the issue.

#. If you are experiencing a problem during the execution of the application, provide
   a script, which allows to quickly reproduce the problem. If the problem is caused by
   a specific set of input data, please, consider providing the corresponding data file
   (preferably, as a TRIQS-compatible :ref:`HDF5-archive <hdf5_tutorial>`).

Thanks!
