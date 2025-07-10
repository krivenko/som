(changelog)=

# Changelog

## [2.1.2](https://github.com/krivenko/som/tree/2.1.2) (2025-07-10)
[Full Changelog](https://github.com/krivenko/som/compare/2.1.1...2.1.2)

- Set up building of Conda packages for released versions.
  The packages are now available in the ['krivenko' channel](
  https://anaconda.org/krivenko/triqs_som).

## [2.1.1](https://github.com/krivenko/som/tree/2.1.1) (2025-06-30)
[Full Changelog](https://github.com/krivenko/som/compare/2.1.0...2.1.1)

- Accept [TRIQS 3.3](https://github.com/TRIQS/triqs/releases/tag/3.3.0) as
  compatible version.
- Add generation of an
  [easyconfig file](https://docs.easybuild.io/writing-easyconfig-files) as
  part of build process.
- Minor coding style improvements and fixes.

## [2.1.0](https://github.com/krivenko/som/tree/2.1.0) (2023-10-03)
[Full Changelog](https://github.com/krivenko/som/compare/2.0.0...2.1.0)

- Port to [TRIQS 3.2](
  https://github.com/TRIQS/triqs/releases/tag/3.2.0).

## [2.0.0](https://github.com/krivenko/som/tree/2.0.0) (2022-06-27)
[Full Changelog](https://github.com/krivenko/som/compare/1.2...2.0.0)

### Major changes and new features

- Complete port to [TRIQS 3.1](
  https://github.com/TRIQS/triqs/releases/tag/3.1.0) and Python 3.
- Implementation of the Stochastic Optimization with Consistent Constraints
  (SOCC) proposed by Goulko et al in [Phys. Rev. B 95, 014102 (2017)](
  https://doi.org/10.1103/PhysRevB.95.014102). It includes three major pieces
  of functionality.
  * The Consistent Constraints update in the Markov chain used to accumulate
    particular solutions;
  * The Consistent Constraints protocol for constructing final solutions out of
    particular solutions;
  * The solution quality assessment technique implemented in a new Python module
    `som.spectral_stats`.
- For consistency with MaxEnt and other stochastic continuation methods, the
  objective function of the optimization problem has been changed to the
  "goodness of fit" $\chi^2$-functional.
- Adoption of the $\chi^2$-functional has made it possible to support
  user-supplied covariance matrices of input data as an alternative to simple
  estimated error bars (credits to @snirgaz for proposing this feature).
- A new family of integral kernels for symmetric fermionic Green's functions has
  been introduced. The corresponding observable is called `FermionGfSymm`.
- The `BosonAutoCorr` kernels have been changed to more closely reproduce
  results of the `BosonCorr` kernels for the same input data. Both kernel
  families are defined on the whole energy axis and expect the same spectrum
  normalization constants now (before one had to divide the constants by 2 for
  `BosonAutoCorr`).
- Projection of an observable onto a real frequency mesh can now be
  performed using binning (enabled by default). In this mode the projected
  observable is integrated over bins centered around points of the mesh.
- Further MPI parallelization.
- Massively reworked online documentation.

### Python API changes

- Following a convention change for TRIQS applications, the Python package of
  SOM has been renamed from `pytriqs.applications.analytical_continuation.som`
  to a laconic `som`.
- Functionality of the `run()` method of `SomCore` has been split among a
  few new methods,
  * `accumulate()` -- accumulate particular solutions;
  * `adjust_f()` -- adjust the number of global updates `F`;
  * `compute_final_solution()` -- construct the final solution using the
    standard SOM protocol;
  * `compute_final_solution_cc()` -- construct the final solution using
    Goulko's SOCC protocol.

  One may still call the deprecated `run()`, which is equivalent to calling
  `accumulate()` + `compute_final_solution()`.
- In recent versions of TRIQS it became impossible to use the `<<` syntax to
  fill GF containers from user-defined Python objects. Furthermore,
  high-frequency tail data was separated from the GF containers. As a result,
  that syntax had to be abandoned in favor of a few free functions.
  * `fill_refreq()` -- fill a real-frequency observable from a computed SOM
    solution;
  * `compute_tail()` -- compute high-frequency tail coefficients from a
    computed SOM solution;
  * `reconstruct()` -- reconstruct input from a computed SOM solution.
- It is now possible to resume accumulation of particular solutions by calling
  `SomCore.accumulate()` multiple times, and to discard all accumulated
  solutions by calling `SomCore.clear()`.
- A handful of new properties and accessor methods have been added to `SomCore`.
- The `rectangle` and `configuration` C++ objects are now exposed as Python
  classes `Rectangle` and `Configuration`. `Configuration` objects can be saved
  to/loaded from HDF5 archives.
- Updated signature of `som.count_good_solutions()` to take both
  `good_chi_abs` and `good_chi_rel` (thresholds on $\chi$ and
  $\chi/\chi_\mathrm{min}$ for a solution to be considered good).
- A new utility function `som.estimate_boson_corr_spectrum_norms()` has been
  added. Given a correlator of boson-like operators $\chi$ defined on any
  supported mesh, it returns a list of spectrum normalization constants
  $\mathcal{N} = \pi \chi(i\Omega = 0)$.

### Build system and developer tools

- Minimum required CMake version has been bumped to 3.12.4.
- Structure of the project has been adjusted to follow conventions established
  by the [app4triqs](https://github.com/TRIQS/app4triqs) application template.
- A `Dockerfile` has been added.
- Files `som.modulefile` and `somvars.sh` are generated and installed as part
  of the build process.
- New benchmarks: `all_kernels`, `binning`, `consistent_constraints`,
  `bosonautocorr` and `fermiongfsymm`.
- The `chi` benchmark has been removed as it depended on the private
  `triqs_ctseg` code.
- Support for C++ static analysis tools `clang-tidy` and `cppcheck` has been
  added.
- A CMake option has been added to link `libsom` and unit tests to Clang
  sanitizers (`AddressSanitizer` and `UndefinedBehaviorSanitizer`).
- C++/Python coding style is enforced with `clang-format` and `flake8`.

## [1.2](https://github.com/krivenko/som/tree/1.2) (2020-03-15)
[Full Changelog](https://github.com/krivenko/som/compare/1.1...1.2)

- Improvements and small fixes in documentation.
- Added a new histogram post-processing function, `count_good_solutions()`.
- Fixed a bug in `update_glue_shift` elementary update.
- New benchmark `all_kernels` and fixes in the `chi` benchmark.
- Added Travis CI config for continuous testing and documentation deployment.
- Minor code improvements.

## [1.1](https://github.com/krivenko/som/tree/1.1) (2017-04-23)
[Full Changelog](https://github.com/krivenko/som/compare/1.0...1.1)

- Massive extension of documentation.
- `adjust_f` mode is disabled by default.
- Fixed a critical bug in `back_transform()`.

## [1.0](https://github.com/krivenko/som/tree/1.0) (2017-03-19)

- First official release.
- Supports 12 integral kernels:
    * 4 observable kinds (`FermionGf`, `BosonCorr`, `BosonAutoCorr`, and
      `ZeroTemp`);
    * 3 input meshes (`imtime`, `imfreq`, `legendre`).
