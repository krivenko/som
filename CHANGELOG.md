# Change Log

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
    * 4 observable kinds (`FermionGf`, `BosonCorr`, `BosonAutoCorr`, and `ZeroTemp`);
    * 3 input meshes (`imtime`, `imfreq`, `legendre`).
