# Cross-language validation against pyfgls

These tools are excluded from the package build (`.Rbuildignore`).

## Golden fixtures (`tests/testthat/fixtures/golden/`)

`make_golden.py` simulates ten scenarios with pyfgls (pinned version in
`MANIFEST.json`), writes the data (`*_data.csv`, shortest round-trip decimals) and
the pyfgls outputs and intermediates (`*_golden.json`). R never re-simulates golden
data. The tolerance policy is in `tests/testthat/helper-golden.R`.

```sh
pip install "pyfgls==0.1.0"
python tools/make_golden.py            # regenerate
python tools/make_golden.py --check    # compare a fresh regeneration with the committed files
```

`.github/workflows/golden.yaml` runs `--check` weekly with the pinned and the latest
pyfgls, then the R golden tests, and opens an issue on drift.

Known pyfgls 0.1.0 issues (fixed in emiliogcs/pyfgls#2, not merged): its Lasso path
(scikit-learn LARS) stops at penalties below float32 epsilon (~1.2e-7), and
coefficients of order 1e-17 can count as selected. The golden CV curves are
therefore compared for penalties of at least 1.2e-7.

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
