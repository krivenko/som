Code architecture changes in SOM 2.0
====================================

``configuration``
-----------------

- New function ``make_nonoverlapping()`` to convert a generic configuration into
  a non-overlapping one.

``solution_worker``
-------------------

- Add a subdirectory ``c++/som/solution_functionals`` containing implementations
  of the following functionals of an elementary solution.

  * ``objective_function`` (change it to follow the definition (4) from [2])
  * ``fit_quality``
  * ``derivative_penalty`` (Eq. (5) from [2])
  * ``amplitude_penalty`` (Eq. (6) from [2])
  * ``model_deviation`` (Eq. (7) from [2])

- Implement the consistent-constraints (CC) updates of [2].

``som_core``
------------

- New public method ``accumulate()`` as a replacement of ``run()``. Can be called
  multiple times to add more basic solutions. Retain the following parameters

  * ``energy_window``
  * ``max_time``
  * ``verbosity``
  * ``t``
  * ``f``
  * ``l``
  * ``random_seed``
  * ``random_name``
  * ``max_rects``
  * ``min_rect_width``
  * ``min_rect_weight``
  * ``distrib_d_max``
  * ``gamma``
  * ``make_histograms``
  * ``hist_max``
  * ``hist_n_bins``

  and add some CC update-specific parameters.

- New method ``clear()`` to remove all accumulated basic solutions and histograms.

- New free function ``adjust_f(SomCore, f_range, params) -> f`` with the following
  parameters

  * ``energy_window``
  * ``max_time``
  * ``verbosity``
  * ``t``
  * ``l``
  * ``kappa``
  * ``random_seed``
  * ``random_name``
  * ``max_rects``
  * ``min_rect_width``
  * ``min_rect_weight``
  * ``distrib_d_max``
  * ``gamma``
  * CC update-specific parameters

- New free function ``adjust_l(SomCore, l_range, params) -> l`` with the following
  parameters

  * ``energy_window``
  * ``max_time``
  * ``verbosity``
  * ``t``
  * ``f``
  * ``good_chi``
  * ``verygood_chi``
  * ``ratio``
  * ``random_seed``
  * ``random_name``
  * ``max_rects``
  * ``min_rect_width``
  * ``min_rect_weight``
  * ``distrib_d_max``
  * ``gamma``
  * ``make_histograms``
  * ``hist_max``
  * ``hist_n_bins``
  * CC update-specific parameters

- ``adjust_f()`` and ``adjust_l()`` should call ``SomCore.clear()`` at the beginning
  of the adjustment procedure but not at the end.

- New method ``compute_final_solution(good_chi)`` (select basis solutions using the
  standard SOM criterion).

- New method ``compute_final_solution_cc(...)`` (selection algorithm described in [1]).

- Add option ``with_binning`` to the ``fill_observable(GfReFreq)`` method.

New module ``som.spectral_stats``
---------------------------------

- New function ``spectral_avg`` to compute ``i_m`` (Eq. (6) of [1]).
  Supports two overloads,

  * ``spectral_avg(SomCore, MeshReFreq, resolution_function)``, regular energy mesh;
  * ``spectral_avg(SomCore, np.array, resolution_function)``, general energy mesh.

- New function ``spectral_disp`` to compute ``\sigma_m^2`` (Eq. (7) of [2]).
  Supports two overloads,

  * ``spectral_disp(SomCore, MeshReFreq, avg, resolution_function)``, regular energy mesh;
  * ``spectral_disp(SomCore, np.array, avg, resolution_function)``, general energy mesh.

- New function ``spectral_corr`` to compute ``\sigma_{m m'}`` (Eq. (8) of [2]).
  Supports two overloads,

  * ``spectral_corr(SomCore, MeshReFreq, avg, resolution_function)``, regular energy mesh;
  * ``spectral_corr(SomCore, np.array, avg, resolution_function)``, general energy mesh.

  ``resolution_function`` is one of ``rectangle``  (Eq. (3) of [1]), ``lorentzian``
   (Eq. (5) of [1]) and ``gaussian``.

References
----------

[1]: O. Goulko et al. Phys. Rev. B 95, 014102 (2017).
[2]: N. V. Prokof'ev and B. V. Svistunov, JETP Lett., 2013, Vol. 97, No. 11, pp. 649-653.
