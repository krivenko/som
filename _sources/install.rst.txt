.. _install:

Installation
============

.. highlight:: bash

Docker
------

A Docker image including the latest version of SOM as well as other TRIQS
applications is available
`here <https://hub.docker.com/repository/docker/ikrivenko/som>`_.

Compiling SOM from source
-------------------------

Prerequisites
*************

The :ref:`TRIQS <triqslibs:welcome>` toolbox version |triqs_branch|
(see :ref:`TRIQS installation instructions <triqslibs:triqs_install>`).
In the following, we will assume that it is installed in the directory
``path_to_triqs``.

Installation steps
******************

#. Download the source code of the latest stable version by cloning the
   ``krivenko/som`` repository from GitHub:

.. substitution-code-block::

   $ git clone -b |som_version| https://github.com/krivenko/som.git som.src

#. Create and move to a new directory where you will compile the code::

   $ mkdir som.build && cd som.build

#. Ensure that your shell contains the TRIQS environment variables by sourcing
   the ``triqsvars.sh`` file from your TRIQS installation::

   $ source path_to_triqs/share/triqs/triqsvars.sh

#. In the build directory call ``cmake``, including any additional custom CMake
   options, see below::

   $ cmake ../som.src

#. Compile the code, run the unit tests and install the application::

   $ make
   $ make test
   $ make install

.. _install_options:

Custom CMake options
********************

The compilation of SOM can be configured using CMake-options:

::

   cmake ../som.src -DOPTION1=value1 -DOPTION2=value2 ...

.. list-table::
    :header-rows: 1
    :widths: 70 30

    * - Options
      - Syntax
    * - Specify an installation path other than path_to_triqs
      - ``-DCMAKE_INSTALL_PREFIX=path_to_som``
    * - Build in Debugging Mode
      - ``-DCMAKE_BUILD_TYPE=Debug``
    * - Enable compilation of shared libraries
      - ``-DBUILD_SHARED_LIBS=ON``
    * - Build a Debian package
      - ``-DBUILD_DEBIAN_PACKAGE=ON``
    * - Disable testing (not recommended)
      - ``-DBuild_Tests=OFF``
    * - Build the documentation locally
      - ``-DBuild_Documentation=ON``
    * - Build without Python support
      - ``-DPythonSupport=OFF``
    * - Initial size of the cache to store computed LHS
      - ``-DCache_Size=0x4000``
    * - Enable extended debugging output (*developers only*)
      - ``-DANALYZE_SOURCES=ON``
    * - Run static analyzer tools (``clang-tidy`` and ``cppcheck``) as part of
        compilation process (*developers only*)
      - ``-DANALYZE_SOURCES=ON``
    * - Compile SOM library with LLVM Address Sanitizer (*developers only*)
      - ``-DASAN=ON``
    * - Compile SOM library with LLVM Undefined Behavior Sanitizer
        (*developers only*)
      - ``-DUBSAN=ON``
