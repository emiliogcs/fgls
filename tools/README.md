# Cross-language validation against pyfgls

These tools are excluded from the package build (`.Rbuildignore`).

## Golden fixtures (`tests/testthat/fixtures/golden/`)

`make_golden.py` simulates ten scenarios with pyfgls (pinned version in
`MANIFEST.json`), writes the data (`*_data.csv`, shortest round-trip decimals) and
the pyfgls outputs and intermediates (`*_golden.json`). R never re-simulates golden
data. The tolerance policy is in `tests/testthat/helper-golden.R`.

```sh
pip install "pyfgls==0.1.1"
python tools/make_golden.py            # regenerate
python tools/make_golden.py --check    # compare a fresh regeneration with the committed files
```

`.github/workflows/golden.yaml` runs `--check` weekly with the pinned and the latest
pyfgls, then the R golden tests, and opens an issue on drift.

`--check` compares every array, and every data column, norm-wise (max |difference| /
max |value|, as the R tests do) at a relative tolerance of 1e-10, not byte for byte:
the results depend on the CPU in the last bits. numpy's SIMD paths for `np.exp`
(AVX-512 or AVX2) change the simulated `y` of s8 and s9 by 1-32 ulps (reproduce with
`NPY_DISABLE_CPU_FEATURES=X86_V4`), and so do OpenBLAS's kernels for the results
(`OPENBLAS_CORETYPE=Haswell`). Iterative solvers (GMA's MA(1) maximum likelihood,
libsvm's SVR) amplify this to ~1e-6, so their outputs are compared at 1e-5.

pyfgls 0.1.1 fixed the Lasso issues of 0.1.0 (path truncated at penalties below
float32 epsilon; coefficients of order 1e-17 counted as selected). One remains in
0.1.1: when the rescaled covariates are tiny, scikit-learn's LARS reports a degenerate
active set and pyfgls falls back to coordinate descent, whose CV errors at very small
penalties are inexact (in the fixtures: s10, psi = 2, penalties below 4e-10, by up to
5e-3). The golden CV curves are therefore still compared for penalties of at least
1.2e-7. A fix is pending in pyfgls; with it, R and pyfgls agree on every CV curve to
~1e-11. Once a pyfgls release contains it, regenerate the fixtures with it
(`make_golden.py`, updating `PYFGLS_PIN`) and remove the mask (`.lars_valid_lambda` in
`tests/testthat/test-golden-hetero.R`).

## Distributional comparison (`compare_python_r.R`, `make_validation.py`)

```sh
Rscript tools/compare_python_r.R simulate 500 4    # data + folds (R RNG) -> validation/cache
python  tools/make_validation.py validation/cache 4
Rscript tools/compare_python_r.R compare 4          # -> validation/cross_language_report.html
                                                    #    inst/extdata/cross_language.rds
```

Both implementations fit the same data sets with the same cross-validation folds.
The summary (`inst/extdata/cross_language.rds`) feeds the validation vignette and
the slow test `test-statistical-crosslang.R`.
