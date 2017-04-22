.. _install:

.. highlight:: bash

Installation
============


Prerequisite
-------------------

#. The :ref:`TRIQS <triqslibs:welcome>` toolbox version 1.4 (see :ref:`TRIQS installation instruction <triqslibs:installation>`).
   In the following, we will suppose that it is installed in the ``path_to_triqs`` directory.

Installation steps
------------------

#. Download the sources of the solver from github::

     $ git clone https://github.com/krivenko/som.git som.src

#. Create an empty build directory where you will compile the code::

     $ mkdir som.build && cd som.build

#. In the build directory call cmake specifying where the TRIQS library is installed::

     $ cmake -DTRIQS_PATH=path_to_triqs ../som.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

.. note:: Be careful with the cmake command above: set TRIQS_PATH, not CMAKE_INSTALL_PREFIX (this variable is only for the TRIQS library)!

.. _install_options:

Customizing your installation: CMake options
--------------------------------------------

+-------------------------------------------------------+---------------------------------+
| Options                                               | Syntax                          |
+=======================================================+=================================+
| Disable testing (not recommended)                     | -DTests=OFF                     |
+-------------------------------------------------------+---------------------------------+
| Build the documentation locally                       | -DBuild_Documentation=ON        |
+-------------------------------------------------------+---------------------------------+
| Initial size of the cache to store computed LHS       | -DCACHE_SIZE=0x4000             |
+-------------------------------------------------------+---------------------------------+
