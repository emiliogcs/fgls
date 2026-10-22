### fgls (R)

**fgls** is the R implementation of the feasible GLS methods of Perron and
González-Coya (2026, *Journal of Econometrics*) for serially correlated errors and
of González-Coya and Perron (2024, *Journal of Econometric Methods*) for
heteroskedasticity of unknown form, separately or together. It mirrors the Python
package [pyfgls](https://pypi.org/project/pyfgls/) and agrees with it to machine
precision.

```r
pak::pak("emiliogcs/fgls")   # install.packages("fgls") once on CRAN
library(fgls)
fit <- fgls(y ~ x, data = df)            # serial correlation + heteroskedasticity
serial_fgls(y ~ x, data = df)            # serial correlation only
lasso_fgls(y ~ x, data = df)             # heteroskedasticity only
```

Documentation: <https://emiliogcs.com/fgls/> · Source:
<https://github.com/emiliogcs/fgls>
