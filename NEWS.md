# SHIP 2.1.0

## Package metadata and release tooling

- Updated package metadata and maintainer contact information.
- Added GitHub project and issue-tracker URLs.
- Added `testthat` edition 3 configuration.
- Added release notes and exclusions for generated check and archive files.

## Robustness and API

- Added validation for input matrices, covariance targets and gene-group lists.
- Added explicit validation of target types in `build.target()`.
- Handled empty, single-variable and zero-correlation cases safely.
- Returned named components from `shrink.estim()`, including numeric `lambda`.

## Documentation and tests

- Modernized examples to avoid `attach()`.
- Regenerated namespace and Rd documentation.
- Added regression tests for target construction and shrinkage estimation.
