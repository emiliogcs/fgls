## Submission

This is a new submission.

`fgls` implements the feasible GLS estimators of Perron and González-Coya (2026,
Journal of Econometrics, <doi:10.1016/j.jeconom.2026.106346>) and González-Coya and
Perron (2024, Journal of Econometric Methods, <doi:10.1515/jem-2023-0007>). It is the
R counterpart of the Python package 'pyfgls' (PyPI).

## Test environments

* local Ubuntu 24.04, R 4.3.3 (reference BLAS): `R CMD check --as-cran` on the built
  tarball: 0 errors, 0 warnings, 1 note (CRAN incoming feasibility: maintainer and URL
  check; the URLs could not be verified from the build machine, which has no access
  to doi.org, orcid.org or github.com).
* GitHub Actions (to be run once the repository is public): ubuntu-latest (R devel,
  release, oldrel-1, 4.1), macOS-latest and windows-latest (R release).
* win-builder (devel, release) and R-hub v2: to be run before submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for the reviewers

* Words flagged as possibly misspelled in DESCRIPTION (González, Coya, Durbin,
  Cochrane, Orcutt) are author names and established method names.
* Examples run in well under 5 seconds each; slow Monte Carlo tests are skipped on
  CRAN (`skip_on_cran()` and an opt-in environment variable). Vignettes use small
  simulated data sets and precomputed simulation results (`inst/extdata`).
* No examples or tests write outside `tempdir()`; the random-number-generator state
  of the session is restored by every function that takes a `seed`.
