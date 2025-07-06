.. _install:

Installation
============

.. highlight:: bash

.. _install_docker:

Docker
------

A Docker image including the latest version of SOM as well as other TRIQS
applications is available
`here <https://hub.docker.com/r/ikrivenko/som/tags>`_.

.. _install_source:

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

.. _install_easybuild:

Compiling SOM from source using EasyBuild
-----------------------------------------

`EasyBuild <https://docs.easybuild.io/>`_ is a software build and installation
framework that allows you to manage (scientific) software on High Performance
Computing (HPC) systems in an efficient way. Please, make sure that the ``eb``
tool is installed on your system by following
`EasyBuild's installation guide <https://docs.easybuild.io/installation/>`_.

SOM 2.x is available starting from EasyBuild version 5.1.1. To install SOM and
its prerequisites (including toolchains, Python, various libraries and the TRIQS
libraries), type::

    $ eb --robot TRIQS-som-2.1.1-foss-2023a.eb

Corresponding environment modules will also be generated, thus a package can be
loaded using::

    $ module load TRIQS-som/2.1.1-foss-2023a-Python-3.11.3

or simply::

    $ module load TRIQS-som

for the most recent version.
