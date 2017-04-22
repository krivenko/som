.. module:: pytriqs.applications.analytical_continuation.som

.. _documentation:

Documentation
=============

.. toctree::
   :maxdepth: 1

   kernels
   parameters
   examples/fermiongf/example
   examples/bosoncorr/example
   examples/bosonautocorr/example
   examples/zerotemp/example

Running SOM to analytically continue input data requires writing a simple Python script.
Details of the script will vary depending on the physical quantity to be continued, and
its representation (function of imaginary time, Matsubara frequencies, or a list of
the Legendre basis coefficients). Nonetheless, a typical script will the following basic parts.

* Import TRIQS and SOM Python modules.

  ::

        # Green's function containers
        from pytriqs.gf.local import *
        # HDFArchive interface to .h5 files
        from pytriqs.archive import HDFArchive
        # HDF5 archives must be modified only by one MPI rank.
        # is_master_node() checks we are on rank 0.
        from pytriqs.utility.mpi import is_master_node

        # Main SOM class
        from pytriqs.applications.analytical_continuation.som import Som

* *Optional:* Load input data from an HDF5 archive.

  ::

        # Open an HDF5 file in read-only mode
        arch = HDFArchive('input.h5', 'h')

        # Read input function
        inp = arch['input']
        # Read importance function S
        S = arch['S']

  This step can be omitted if the input data comes from a different source.
  For instance, it could be loaded from text files or generated in the very
  same script by a quantum Monte-Carlo impurity solver.

* Construct `Som` object. Here you must specify the kind of the physical
  observable in question, or equivalently, the integral kernel in the analytical
  continuation equation.

  ::

        # Create Som object
        cont = Som(inp, S, kind = "FermionGf", norms = norms)

  `inp` and `S` must be Green's function containers of the same type (one of
  `GfImTime`, `GfImFreq` or `GfLegendre`), defined on the same mesh and having
  the same target shape.

  Currently supported observable kinds are `FermionGf` (fermionic Green's function
  or self-energy), `BosonCorr` (correlator of a boson-like operator with its Hermitian
  conjugate), `BosonAutoCorr` (correlator of a Hermitian operator with itself), and
  `ZeroTemp` (correlator computed in Matsubara formalism but at zero temperature).
  Refer to ':ref:`kernels`' for a detailed description of the implemented observable kinds.

  If target shape of `inp` is bigger than 1x1, SOM will only construct analytical
  continuation of the diagonal matrix elements.

  `norms` is an optional one-dimensional NumPy array containing expected norms of the
  spectral functions to be found, one positive real number per one diagonal element of `inp`.
  Instead of setting all elements of `norms` to the same constant `x`, one may simply
  pass `norms = x`.
  By default all norms are set to 1, which is correct for the Green's functions of
  the fermions. However, adjustments would normally be needed for self-energies and
  bosonic correlators/autocorelators.

* Set simulation parameters.

  ::

        # Dictionary contaning all run() parameters
        run_params = {}

        # SOM will construct a spectral function within this energy window
        run_params['energy_window'] = (-5.0,5.0)

        # Number of particular solutions to calculate. It controls smoothness of
        # the final solution, and should in practice be chosen from 100-10000 range
        # (computation time scales linearly with it).
        run_params['l'] = 1000

  These two parameters are required for any simulation. Other useful parameters
  include `verbosity` and setting of the underlying Markov chain algorithm. You
  can find a detailed discussion of the latter in Sections 3 and 4 of
  :download:`Mishchenko's lecture <notes/som.pdf>`.

  ::

        # Set verbosity level to the maximum on MPI rank 0, and silence all other ranks
        run_params['verbosity'] = 3 if is_master_node() else 0

        # Do not auto-adjust the number of particular solutions to accumulate
        # (default and recommended)
        run_params['adjust_l'] = False

        # Do not auto-adjust the number of global updates (default and recommended)
        run_params['adjust_f'] = False

        # Number of global updates per particular solution
        run_params['f'] = 200

        # Number of local updates per global update
        run_params['t'] = 500

        # Accumulate histogram of the objective function values,
        # useful to analyse quality of solution
        run_params['make_histograms'] = True


  See ':ref:`parameters`' for a full list of accepted parameters.

* Run simulation.

  ::

        cont.run(**run_params)

* Extract solution and reconstructed input.

  SOM internally represents spectral functions
  as sums of rectangles. This representation is accessible as `cont.solutions`. In most
  cases one is interested in the resulting observable defined on a real frequency mesh.

  ::

        sol = cont.solutions
        # 'sol' is now a list of spectral functions, len(sol) == len(inp.indices)
        # Each element of 'sol' represents a sum of rectangles; its type is
        # list((rect_center,rect_width,rect_height))

        # Evaluate the solution on a frequency mesh.
        # The tail coefficients will be written into f_w too.
        #
        # NB: we can use *any* energy window at this point, not necessarily that from run_params
        f_w = GfReFreq(window = (-5.0, 5.0), n_points = 1000, indices = inp.indices)
        f_w << cont

        # Imaginary time/frequency/Legendre data reconstructed from the solution
        rec = inp.copy()
        rec << cont

  It is neccessary (but not enough) to have `rec` close to `inp` to ensure correctness
  of the solution.

* *Optional:* Save results to an HDF5 archive.

  ::

        # On master node, save results to an archive
        if mpi.is_master_node():
            with HDFArchive("results.h5",'w') as ar:
                ar['inp'] = inp
                ar['f_w'] = f_w
                ar['rec'] = rec
                # Save histograms for further analysis
                ar['histograms'] = cont.histograms

There are a few :ref:`examples <documentation>` explaining how to run SOM for specific observable kinds.

