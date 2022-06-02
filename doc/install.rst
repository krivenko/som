.. _install:

.. highlight:: bash

Installation
============

Prerequisite
-------------------

#. The :ref:`TRIQS <triqslibs:welcome>` toolbox version 2.2.x (see :ref:`TRIQS installation instruction <triqslibs:installation>`).
   In the following, we will assume that it is installed in the ``path_to_triqs`` directory.

Installation steps
------------------

#. Download the sources of the solver from github::

    $ git clone https://github.com/krivenko/som.git som.src

#. Create an empty build directory where you will compile the code::

    $ mkdir som.build && cd som.build

#. Ensure that your shell contains the TRIQS environment variables by sourcing the ``triqsvars.sh`` file from your TRIQS installation::

    $ source path_to_triqs/share/triqsvarsh.sh

#. In the build directory call ``cmake``, including any additional custom CMake options, see below::

    $ cmake ../som.src

#. Compile the code, run the tests and install the application::

    $ make
    $ make test
    $ make install

.. _install_options:

Customizing your installation: CMake options
--------------------------------------------

The compilation of SOM can be configured using CMake-options::

    cmake ../som.src -DOPTION1=value1 -DOPTION2=value2 ... ../som.src

+-------------------------------------------------------------+------------------------------------+
| Options                                                     | Syntax                             |
+=============================================================+====================================+
| Specify an installation path other than path_to_triqs       | -DCMAKE_INSTALL_PREFIX=path_to_som |
+-------------------------------------------------------------+------------------------------------+
| Disable testing (not recommended)                           | -DBuild_Tests=OFF                  |
+-------------------------------------------------------------+------------------------------------+
| Build the documentation locally                             | -DBuild_Documentation=ON           |
+-------------------------------------------------------------+------------------------------------+
| Initial size of the cache to store computed LHS             | -DCache_Size=0x4000                |
+-------------------------------------------------------------+------------------------------------+
| Build in Debugging Mode                                     | -DCMAKE_BUILD_TYPE=Debug           |
+-------------------------------------------------------------+------------------------------------+
| Enable extended debugging output (*developers only*)        | -DExt_Debug=ON                     |
+-------------------------------------------------------------+------------------------------------+
| Run static analyzer tools (``clang-tidy`` and ``cppcheck``) | -DANALYZE_SOURCES=ON               |
| as part of compilation process (*developers only*)          |                                    |
+-------------------------------------------------------------+------------------------------------+
| Compile SOM library with LLVM Address Sanitizer             | -DASAN=ON                          |
| (*developers only*)                                         |                                    |
+-------------------------------------------------------------+------------------------------------+
| Compile SOM library with LLVM Undefined Behavior Sanitizer  | -DUBSAN=ON                         |
| (*developers only*)                                         |                                    |
+-------------------------------------------------------------+------------------------------------+
