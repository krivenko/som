.. _tutorial:

Tutorial
========

.. currentmodule:: som

Running SOM to analytically continue input data requires writing a simple Python
script. Details of the script will vary depending on the physical quantity to be
continued, and its representation (function of imaginary time, Matsubara
frequencies, or a list of Legendre basis coefficients). Nonetheless, a typical
script will have the following basic parts.

* Import TRIQS and SOM Python modules.

  ::

        # Green's function containers
        from triqs.gf import *
        # HDFArchive interface to .h5 files
        from h5 import HDFArchive
        # HDF5 archives must be modified only by one MPI rank.
        # is_master_node() checks whether we are on rank 0.
        from triqs.utility.mpi import is_master_node

        # Main SOM class
        from som import Som

* *Optional:* Load input data from an
  :ref:`HDF5 archive <triqslibs:hdf5_tutorial>`.

  ::

        # Open an HDF5 file in read-only mode
        arch = HDFArchive('input.h5', 'r')

        # Read input function
        inp = arch['input']

        # Read estimated error bars or a full covariance matrix of
        # the input data
        error_bars = arch['error_bars']
        cov_matrices = arch['cov_matrices']

  This step can be omitted if the input data comes from a different source.
  For instance, it could be loaded from text files or generated in the very
  same script by a quantum Monte-Carlo impurity solver.

* Construct the :class:`Som` object. Here, you must specify the kind of the
  physical observable in question.

  ::

        # Create Som object using the estimated error bars
        cont = Som(inp, error_bars, kind = "FermionGf", norms = norms)

  ::

        # Create Som object using the covariance matrices
        cont = Som(inp,
                   cov_matrices,
                   kind = "FermionGf",
                   norms = norms,
                   filtering_levels = 1e-5)

  ``inp`` and ``error_bars`` must be Green's function containers of
  the same type :class:`triqs.gf.gf.Gf`, defined on the same mesh
  (:class:`triqs.gf.meshes.MeshImTime`, :class:`triqs.gf.meshes.MeshImFreq` or
  :class:`triqs.gf.meshes.MeshLegendre`) and having the same square target
  shape.

  Currently supported observable kinds are

  * :ref:`FermionGf <fermiongf>` - fermionic Green's function or self-energy;
  * :ref:`FermionGfSymm <fermiongfsymm>` - fermionic Green's function or
    self-energy with enforced particle-hole symmetry.
  * :ref:`BosonCorr <bosoncorr>` - bosonic Green's function, transversal
    magnetic susceptibility, general dynamical correlator of the form
    :math:`\chi_{OO^\dagger}(\tau) = \langle \mathbf{T}_\tau \hat O(\tau)
    \hat O^\dagger(0)\rangle`, where :math:`\hat O` is a boson-like
    (fermion-number-conserving) operator.
  * :ref:`BosonAutoCorr <bosonautocorr>` - charge susceptibility,
    longitudinal magnetic susceptibility, optical conductivity,
    general dynamical correlator of the form
    :math:`\chi_{OO}(\tau) = \langle \mathbf{T}_\tau \hat O(\tau)
    \hat O(0)\rangle`.
  * :ref:`ZeroTemp <zerotemp>` - dynamical correlator computed in the
    Matsubara formalism but at zero temperature.

  Refer to ':ref:`observables`' for a detailed description of the implemented
  observable kinds.

  If target shape of ``inp`` is :math:`M{\times}M` with :math:`M>1`,
  SOM will solve an independent analytic continuation problem for each of its
  diagonal matrix elements.

  ``norms`` is an optional one-dimensional NumPy array containing :math:`M`
  expected  norms of the spectral functions to be found, one positive real
  number per one diagonal element of ``inp``. Instead of setting all elements of
  ``norms`` to the same constant ``x``, one may simply pass ``norms = x``.
  By default, all norms are set to 1.0 for the ``FermionGf``, ``FermionGfSymm``
  and ``ZeroTemp`` observables, which is the correct choice for the fermionic
  Green's functions. However, adjustments would normally be needed for the
  fermionic self-energies. ``BosonCorr`` and ``BosonAutoCorr`` always require
  the norms to be explicitly specified. If those normalization constants are
  not known from some analytical considerations, they can be estimated from
  the input data itself using :func:`estimate_boson_corr_spectrum_norms()`.

  The :ref:`covariance matrices <cov_matrix>` :math:`\hat\Sigma` must be
  provided as a container of type
  ``Gf(mesh=MeshProduct(inp.mesh, inp.mesh), target_shape=[M])``.
  Each target element of this container holds one covariance matrix. It is
  also advised to specify positive :ref:`filtering levels <cov_matrix_filtered>`
  via the ``filtering_levels`` argument.

* Set parameters for accumulation of
  :ref:`particular solutions <particular_solutions>`.

  ::

        # Dictionary containing all accumulate() parameters
        acc_params = {}

        # SOM will construct a spectral function within this energy window
        acc_params['energy_window'] = (-5.0, 5.0)

        # Number of particular solutions to accumulate. It controls smoothness
        # of the final solution, and should in practice be chosen from the
        # 100-10000 range (computation time scales linearly with it).
        acc_params['l'] = 1000

  These two parameters are required for any simulation. Other useful parameters
  include ``verbosity`` and settings of the underlying Markov chain algorithm.
  You can find a detailed discussion of the latter in Sections 3 and 4 of
  `A. S. Mishchenko's lecture
  <https://www.cond-mat.de/events/correl12/manuscripts/mishchenko.pdf>`_.

  ::

        # Set verbosity level to the maximum on MPI rank 0, and silence all
        # other ranks.
        acc_params['verbosity'] = 3 if is_master_node() else 0

        # Do not auto-adjust the number of particular solutions to accumulate
        # (default and recommended).
        acc_params['adjust_l'] = False

        # Do not auto-adjust the number of global updates (default and
        # recommended).
        acc_params['adjust_f'] = False

        # Number of global updates per particular solution.
        # Bigger values improve ergodicity of the Markov chain algorithm
        # (computation time scales linearly with 'f').
        acc_params['f'] = 200

        # Number of local updates per global update.
        # Has a similar effect on the algorithm as 'f', and scales linearly too.
        acc_params['t'] = 500

        # Enable the Consistent Constraints updates to speed-up accumulation.
        acc_params['cc_update'] = True

        # Accumulate histogram of the objective function values \chi,
        # useful to analyse quality of solution.
        acc_params['make_histograms'] = True


  See the documentation of the method :class:`Som.accumulate()` for a full list
  of accepted parameters.

* Accumulate particular solutions.

  ::

        cont.accumulate(**acc_params)

  Normally, this is the most time consuming step. Calling
  :class:`Som.accumulate()` multiple times will incrementally extend the pool
  of the accumulated solutions.

* Construct the :ref:`final solution <final_solution>` out of the accumulated
  particular solutions. One can choose to use either the standard SOM procedure
  and average selected particular solutions with equal weights,

  ::

      cont.compute_final_solution()

  or the more advanced :ref:`SOCC procedure <final_solution_cc>` introduced in
  SOM 2.0.

  ::

      cont.compute_final_solution_cc()

  Behavior of :func:`Som.compute_final_solution_cc()` can be fine-tuned via
  the many keyword arguments it accepts.

* Extract the final solution,
  :ref:`recover the real-frequency version <recovery>` of ``inp``,
  its :ref:`tail expansion coefficients <triqslibs:gf_tail>` and
  reconstruct the input.

  SOM internally represents spectral functions (configurations) as sums of
  rectangles. The final solutions are accessible as a list of
  :class:`Configuration` objects via an attribute :data:`Som.solutions`.

  ::

        sol = cont.solutions
        # 'sol' is now a list of spectral functions,
        # len(sol) == len(inp.indices)

        # Recover the real-frequency version of 'inp'
        #
        # NB: we can use *any* energy window at this point, not necessarily that
        # from 'acc_params'
        w_mesh = MeshReFreq(window=(-5.0, 5.0), n_w=1000)
        f_w = Gf(mesh=w_mesh, indices=inp.indices)
        fill_refreq(f_w, cont)

        # High frequency expansion (tail) coefficients up to order 'max_order'
        # can be computed by function 'compute_tail()' and returned as a
        # 3-dimensional complex NumPy array.
        # The first index of the array is the zero-based expansion order.
        tail = compute_tail(5, cont) # max_order = 5

        # Imaginary time/frequency/Legendre data reconstructed from the solution
        rec = inp.copy()
        reconstruct(rec, cont)

  It is necessary (but not sufficient) to have ``rec`` close to ``inp``
  to ensure good quality of the final solution.

* *Optional:* Save results to an HDF5 archive.

  ::

        # On master node, save results to an archive
        if mpi.is_master_node():
            with HDFArchive("results.h5", 'w') as ar:
                ar['inp'] = inp
                ar['sol'] = sol
                ar['f_w'] = f_w
                ar['tail'] = tail
                ar['rec'] = rec
                # Save histograms for further analysis
                ar['histograms'] = cont.histograms

* *Optional*: Study :ref:`statistical properties <spectral_stats>` of the
  accumulated particular solutions.

  ::

      from som.spectral_stats import spectral_avg, spectral_disp, spectral_corr

      # Compute statistical characteristics of the accumulated solution ensemble
      # for the 0-th diagonal matrix element of 'inp' on the real-frequency
      # mesh 'w_mesh' and using the rectangular resolution function.

      # Averages over the ensemble
      avg = spectral_avg(cont, 0, f_w.mesh, "rectangle")
      # Dispersions of the ensemble
      disp = spectral_disp(cont, 0, f_w.mesh, avg, "rectangle")
      # Two-point correlations of the ensemble
      corr = spectral_corr(cont, 0, f_w.mesh, avg, "rectangle")

There are a few :ref:`examples <examples>` explaining how to run SOM for
various observables and how to use its more advanced features.
