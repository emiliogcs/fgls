#!/usr/bin/env python3
"""Generate the golden fixtures used by the R test-suite.

Each scenario writes the data it was computed on (CSV, shortest round-trip decimal
representation of every float64, so R reads the *same* numbers) and the pyfgls
outputs and intermediates (JSON). R never re-simulates golden data.

Usage (from the package root, with the pinned pyfgls installed):

    python tools/make_golden.py [--out tests/testthat/fixtures/golden]

The pyfgls version, the numerical library versions and the seeds are recorded in
``MANIFEST.json``. Fixture sizes are kept small (< 500 KB in total).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np

import pyfgls
from pyfgls._utils import ar_roots
from pyfgls.serial.fgls import estimate_serial
from pyfgls.serial.lagselect import resolve_kmax
from pyfgls.serial.ma import estimate_gma
from pyfgls.serial.variance import serial_variance
from pyfgls.simulate import combined_dgp, hetero_dgp, serial_dgp

PYFGLS_PIN = "0.1.1"
SCENARIOS: dict[str, Any] = {}


def scenario(fn):  # registry decorator
    SCENARIOS[fn.__name__] = fn
    return fn


# ----------------------------------------------------------------------------- output
def _f(v: float) -> str:
    return repr(float(v))


def write_csv(path: Path, cols: dict[str, np.ndarray]) -> None:
    names = list(cols)
    n = len(cols[names[0]])
    with path.open("w", newline="\n") as fh:
        fh.write(",".join(names) + "\n")
        for i in range(n):
            fh.write(",".join(_f(cols[c][i]) for c in names) + "\n")


def jsonable(x: Any) -> Any:
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        if np.iscomplexobj(x):
            return {"re": x.real.tolist(), "im": x.imag.tolist()}
        return x.tolist()
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


# ------------------------------------------------------------------- serial helpers
def serial_block(y, X, names, *, ar_order=None, kmax=None, lagged_dep_idx=(), wallis="auto"):
    """Every intermediate of the PG26 estimator, for the robust and the other variances."""
    T, q = X.shape
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        est = estimate_serial(
            y,
            X,
            names,
            has_const=True,
            ar_order=ar_order,
            kmax=kmax,
            lagged_dep_idx=lagged_dep_idx,
            wallis_regressor=wallis,
            cov_type="robust",
        )
    fr = est.filtered
    var = est.variance
    J = None if est.wallis is None else est.wallis.J
    other = {}
    for ct in ("naive", "HC0", "HC1", "HC2", "HC3"):
        other[ct] = serial_variance(fr.Xf, fr.ehat, fr.uhat, est.design, J=J, cov_type=ct).cov
    sel = est.selection
    out: dict[str, Any] = {
        "kstar": est.order,
        "kmax": None if sel is None else sel.kmax,
        "selected_by_ic": None if sel is None else sel.selected_by_ic,
        "bic_path": None if sel is None else sel.values,
        "durbin_names": list(est.design.names),
        "dropped_kstar": list(est.design.dropped),
        "beta_durbin": est.durbin.beta,
        "durbin_coef": est.durbin.coef,
        "durbin_bse": est.durbin.bse,
        "rho_hat": est.rho,
        "ar_roots_mod": np.sort(np.abs(ar_roots(est.rho))),
        "beta_fgls": fr.beta,
        "e_hat": fr.ehat,
        "u_hat": fr.uhat,
        "n_eff": fr.n_eff,
        "Ghat": var.G,
        "Sigma_psi": var.Sigma,
        "Qhat": fr.Xf.T @ fr.Xf / fr.n_eff,
        "sigma2_e": var.sigma2,
        "AVar_naive": var.naive_avar,
        "AVar_R": var.avar,
        "vcov": var.cov,
        "vcov_other": other,
        "nonexogeneity_by_coef": var.nonexogeneity_by_coef,
        "nonexogeneity_share": var.nonexogeneity_share,
    }
    if sel is not None:
        # dropped columns of the common-sample design used by the lag selection
        from pyfgls.serial.durbin import durbin_design

        out["dropped_kmax"] = list(durbin_design(y, X, names, sel.kmax, start=sel.kmax).dropped)
    if est.wallis is not None:
        out["wallis"] = {
            "rho": est.wallis.rho,
            "J": est.wallis.J,
            "used": list(est.wallis.used),
            "ratios": {k: v for k, v in est.wallis.ratios.items()},
        }
    return out


SLIM_KEYS = (
    "method", "beta_fgls", "weights", "vcov_eq6_8", "vcov_types", "df_hat", "fitted_logvar",
    "nonzero_set", "psi_star", "lambda_star", "intercept", "gamma_hat", "svr", "rho_hat",
    "vcov", "vcov_other", "kstar", "durbin_coef", "beta_durbin", "nonexogeneity_share",
)


def slim(block: dict[str, Any]) -> dict[str, Any]:
    """Secondary results: keep only the essential outputs (fixture size budget)."""
    return {k: v for k, v in block.items() if k in SLIM_KEYS}


def public_api_check(y, X, names, **kw) -> dict[str, Any]:
    """End-to-end pyfgls result through the public API (sanity check of the internals)."""
    import pandas as pd

    df = pd.DataFrame(X, columns=names)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = pyfgls.SerialFGLS(pd.Series(y, name="y"), df, **kw).fit()
    return {"params": res.params.to_numpy(), "bse": res.bse.to_numpy(), "ar_order": res.ar_order}


# ------------------------------------------------------------------- hetero helpers
def hetero_block(y, X, Z, z_names, seed, *, method="lasso", cv_folds=10, lambda_rule="min"):
    """Every intermediate of the GP24 estimator (adaptive Lasso by default)."""
    from pyfgls._linalg import hat_diag
    from pyfgls.hetero.adaptive_lasso import make_folds, ridge_cv, standardize
    from pyfgls.hetero.fgls import estimate_hetero
    from pyfgls.hetero.variance import hetero_covariance, weighted_leverage
    from pyfgls.hetero.weights import log_squared_residuals

    T = len(y)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        est = estimate_hetero(
            y, X, Z, z_names, method=method, random_state=seed, cv_folds=cv_folds,
            lambda_rule=lambda_rule,
        )
    u0 = est.u_unweighted
    r = log_squared_residuals(u0, est.delta)
    label = method if est.fit is not None else ("ols" if method == "none" else "known")
    covs = {
        ct: hetero_covariance(
            ct, X=X, v=est.v, u_fgls=est.u_weighted, u_ols=u0, h_ols=est.h_unweighted,
            df=est.df, method=label,
        )
        for ct in ("robust", "nonrobust", "HC0", "HC1", "HC2", "HC3")
    }
    out: dict[str, Any] = {
        "method": method,
        "beta_ols": est.beta_unweighted,
        "u_ols": u0,
        "h_ols": est.h_unweighted,
        "response": r,
        "weights": est.v,
        "h_fgls": weighted_leverage(X, est.v),
        "beta_fgls": est.beta,
        "u_fgls": est.u_weighted,
        "df_hat": est.df,
        "vcov_eq6_8": est.cov,
        "vcov_types": covs,
    }
    fit = est.fit
    if fit is None:
        return out
    out.update(
        {
            "fitted_logvar": fit.fitted,
            "intercept": fit.intercept,
            "gamma_hat": fit.coef,
            "nonzero_set": list(fit.selected),
            "z_names": list(z_names),
        }
    )
    tun = fit.tuning
    if method == "lasso":
        folds = make_folds(T, min(cv_folds, T), "random", seed)
        Zs, mu, sd, keep = standardize(Z)
        g_ridge, lam_ridge = ridge_cv(Zs, r, folds)
        psis = sorted(tun["cv_path"])
        psi, lam = tun["psi"], tun["lambda"]
        a = np.abs(g_ridge)
        g_std = fit.coef[keep] * sd[keep]
        with np.errstate(divide="ignore"):
            w = np.where(a > 0, a ** (-psi), np.inf)
        pen = float(np.sum(np.where(g_std != 0, w * np.abs(g_std), 0.0)))
        rss = float(np.sum((r - fit.fitted) ** 2))
        out.update(
            {
                "foldid": folds + 1,
                "ridge_coef": g_ridge,
                "lambda_ridge": lam_ridge,
                "keep": keep + 1,
                "z_mean": mu,
                "z_sd": sd,
                "psi_grid": psis,
                # each grid is lam_max * 10^seq(0, log10(ratio), length = n_lambda)
                "cv_lambda_range": [
                    [tun["cv_path"][p][0][0], tun["cv_path"][p][0][-1]] for p in psis
                ],
                "n_lambda": len(tun["cv_path"][psis[0]][0]),
                "cv_error_matrix": [tun["cv_path"][p][1] for p in psis],
                "psi_star": psi,
                "lambda_star": lam,
                "lambda_paper": T * lam,
                "cv_mse": tun["cv_mse"],
                "objective": 0.5 * rss + T * lam * pen,
            }
        )
    elif method == "svr":
        out.update({"svr": {k: tun[k] for k in ("C", "epsilon", "cv_mse", "n_support")}})
        out["foldid"] = make_folds(T, min(cv_folds, T), "random", seed) + 1
    return out


def lasso_public(y, X, x_names, Z, z_names, seed) -> dict[str, Any]:
    import pandas as pd

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = pyfgls.LassoFGLS(
            pd.Series(y, name="y"), pd.DataFrame(X, columns=x_names),
            z=pd.DataFrame(Z, columns=z_names), random_state=seed,
        ).fit()
    return {"params": res.params.to_numpy(), "bse": res.bse.to_numpy(),
            "selected": list(res.selected_z)}


def combined_block(y, X, x_names, seed, *, first_stage_correction=True, cv_folds=10):
    """F-C / F-CH ([PG26] Supplement S-4) with z = "auto" built from the raw regressors."""
    from pyfgls.combined import estimate_combined
    from pyfgls.hetero.adaptive_lasso import make_folds
    from pyfgls.hetero.features import make_dictionary
    from pyfgls.hetero.variance import weighted_leverage

    Z, zn = make_dictionary(X, x_names)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        est = estimate_combined(
            y, X, x_names, Z, zn, random_state=seed, cv_folds=cv_folds,
            first_stage_correction=first_stage_correction,
        )
    ser, het = est.serial, est.hetero
    fr = ser.filtered
    out: dict[str, Any] = {
        "kstar": ser.order,
        "rho_hat": ser.rho,
        "beta_FC": est.beta_fc,
        "beta_FCH": est.beta_fch,
        "vcov_FC": est.cov_fc,
        "vcov_FCH": est.cov_fch,
        "weights": het.v,
        "df_hat": het.df,
        "h_fc": weighted_leverage(fr.Xf, np.ones(fr.n_eff)),
        "h_fch": weighted_leverage(fr.Xf, het.v),
        "z_names": list(zn),
        "first_stage_correction": first_stage_correction,
    }
    if het.fit is not None:
        tun = het.fit.tuning
        out.update(
            {
                "nonzero_set": list(het.fit.selected),
                "psi_star": tun.get("psi"),
                "lambda_star": tun.get("lambda"),
                "fitted_logvar": het.fit.fitted,
                "foldid": make_folds(fr.n_eff, min(cv_folds, fr.n_eff), "random", seed) + 1,
            }
        )
    return out


def benchmark_block(y, X, names) -> dict[str, Any]:
    """pyfgls.benchmarks on the same data: OLS (HAC QS, NW, HC3), CO, Durbin, infeasible GLS."""
    import pandas as pd

    from pyfgls import benchmarks as bm

    Xd = pd.DataFrame(X, columns=names)
    ys = pd.Series(y, name="y")
    out: dict[str, Any] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for ct in ("HAC", "NW", "HC3"):
            r = bm.ols(ys, Xd, cov_type=ct)
            out[f"ols_{ct}"] = {"params": r.params.to_numpy(), "cov": r.cov.to_numpy(),
                                "label": r.cov_type}
        co = bm.cochrane_orcutt(ys, Xd)
        out["co"] = {"params": co.params.to_numpy(), "cov": co.cov.to_numpy(),
                     "rho": co.rho.to_numpy()}
        coi = bm.cochrane_orcutt(ys, Xd, iterate=True)
        out["co_iter"] = {"params": coi.params.to_numpy(), "rho": coi.rho.to_numpy()}
        du = bm.durbin(ys, Xd)
        out["durbin"] = {"params": du.params.to_numpy(), "cov": du.cov.to_numpy(),
                         "rho": du.rho.to_numpy()}
        ig = bm.infeasible_gls(ys, Xd, rho=[0.8])
        out["gls_rho"] = {"params": ig.params.to_numpy(), "cov": ig.cov.to_numpy()}
    return out


def diagnostics_block(y, X, names) -> dict[str, Any]:
    """pyfgls.diagnostics: pre-determinedness LM test (both scales) and ADF statistics."""
    import pandas as pd

    from pyfgls import diagnostics as dg

    Xd = pd.DataFrame(X, columns=names)
    out: dict[str, Any] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for sc in ("innovations", "errors"):
            t = dg.test_predetermined(pd.Series(y, name="y"), Xd, max_lag=4, scale=sc)
            out[f"predetermined_{sc}"] = {
                "LM": t.table["LM"].to_numpy(), "df": t.table["df"].to_numpy(),
                "pvalue": t.table["p-value"].to_numpy(), "ar_order": t.ar_order,
            }
        a = dg.adf_test(y)
        out["adf_y"] = {"statistic": a.statistic, "lags": a.lags, "nobs": a.nobs}
        a2 = dg.adf_test(X[:, 1])
        out["adf_x"] = {"statistic": a2.statistic, "lags": a2.lags, "nobs": a2.nobs}
    return out


def gp24_z(x):
    """The GP24 Section 5 candidate set (x, [log x]^2, x^2, cos x, cos 2x)."""
    names = ["x", "log(x)^2", "x^2", "cos(x)", "cos(2x)"]
    return np.column_stack([x, np.log(x) ** 2, x**2, np.cos(x), np.cos(2.0 * x)]), names


# -------------------------------------------------------------------------- scenarios
@scenario
def s1_ar1_nonexog(out: Path) -> dict[str, Any]:
    """AR(1) errors (0.8), non-exogenous regressor (gamma = 0.5), PG26 Table-1 design."""
    sim = serial_dgp(500, ar=(0.8,), rho_x=0.0, gamma=0.5, random_state=20260101)
    names = ["(Intercept)", "x"]
    write_csv(out / "s1_data.csv", {"y": sim.y, "x": sim.X[:, 1]})
    res = serial_block(sim.y, sim.X, names)
    res["fixed_order_3"] = slim(serial_block(sim.y, sim.X, names, ar_order=3))
    res["public"] = public_api_check(sim.y, sim.X, names)
    res["benchmarks"] = benchmark_block(sim.y, sim.X, names)
    res["diagnostics"] = diagnostics_block(sim.y, sim.X, names)
    return {"formula": "y ~ x", "T": 500, "seed": 20260101, "results": res}


@scenario
def s2_ar2_exog(out: Path) -> dict[str, Any]:
    """AR(2) errors (1.34, -0.42), exogenous persistent regressor."""
    sim = serial_dgp(200, ar=(1.34, -0.42), rho_x=0.8, gamma=0.0, burn=100, random_state=20260102)
    names = ["(Intercept)", "x"]
    write_csv(out / "s2_data.csv", {"y": sim.y, "x": sim.X[:, 1]})
    res = serial_block(sim.y, sim.X, names)
    res["public"] = public_api_check(sim.y, sim.X, names)
    return {"formula": "y ~ x", "T": 200, "seed": 20260102, "results": res}


@scenario
def s3_ma1(out: Path) -> dict[str, Any]:
    """MA(1) errors (theta = -0.4): AR-approximation FGLS and GMA."""
    sim = serial_dgp(200, ar=(), ma=(-0.4,), rho_x=0.8, gamma=0.0, burn=100, random_state=20260103)
    names = ["(Intercept)", "x"]
    write_csv(out / "s3_data.csv", {"y": sim.y, "x": sim.X[:, 1]})
    res = serial_block(sim.y, sim.X, names)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = estimate_gma(sim.y, sim.X, names)
    res["gma"] = {
        "theta": g.theta,
        "durbin_order": g.durbin_order,
        "beta": g.beta,
        "vcov": g.cov,
        "e_hat": g.ehat,
    }
    return {"formula": "y ~ x", "T": 200, "seed": 20260103, "results": res}


@scenario
def s4_white_noise(out: Path) -> dict[str, Any]:
    """White-noise errors: the BIC should pick k* = 0 (FGLS = OLS)."""
    sim = serial_dgp(200, ar=(), ma=(), rho_x=0.5, gamma=0.0, burn=100, random_state=20260104)
    names = ["(Intercept)", "x"]
    write_csv(out / "s4_data.csv", {"y": sim.y, "x": sim.X[:, 1]})
    res = serial_block(sim.y, sim.X, names)
    return {"formula": "y ~ x", "T": 200, "seed": 20260104, "results": res}


@scenario
def s5_deterministic(out: Path) -> dict[str, Any]:
    """Constant + linear trend + quarterly dummy: collinear lagged block in the Durbin regression."""
    rng = np.random.default_rng(20260105)
    T, burn = 200, 100
    e = rng.standard_normal(T + burn)
    v = rng.standard_normal(T + burn)
    u = np.zeros(T + burn)
    x = np.zeros(T + burn)
    for t in range(1, T + burn):
        u[t] = 0.6 * u[t - 1] + e[t]
        x[t] = 0.5 * x[t - 1] + v[t]
    u, x = u[burn:], x[burn:]
    trend = np.arange(1.0, T + 1.0)
    dum = (np.arange(T) % 4 == 0).astype(float)
    y = 1.0 + 0.5 * x + 0.02 * trend + 0.5 * dum + u
    X = np.column_stack([np.ones(T), x, trend, dum])
    names = ["(Intercept)", "x", "trend", "dum"]
    write_csv(out / "s5_data.csv", {"y": y, "x": x, "trend": trend, "dum": dum})
    res = serial_block(y, X, names)
    return {"formula": "y ~ x + trend + dum", "T": T, "seed": 20260105, "results": res}


@scenario
def s6_lagged_dep(out: Path) -> dict[str, Any]:
    """Lagged dependent variable as regressor: Wallis recovery of rho (PG26 Section 4.2)."""
    rng = np.random.default_rng(20260106)
    T, burn = 400, 100
    n = T + burn + 1
    e = rng.standard_normal(n)
    v = rng.standard_normal(n)
    u = np.zeros(n)
    x = np.zeros(n)
    yy = np.zeros(n)
    for t in range(1, n):
        u[t] = 0.5 * u[t - 1] + e[t]
        x[t] = 0.6 * x[t - 1] + v[t]
        yy[t] = 0.5 + 0.4 * yy[t - 1] + 1.0 * x[t] + u[t]
    y = yy[burn + 1 :]
    ylag = yy[burn:-1]
    xx = x[burn + 1 :]
    X = np.column_stack([np.ones(T), ylag, xx])
    names = ["(Intercept)", "ylag1", "x"]
    write_csv(out / "s6_data.csv", {"y": y, "ylag1": ylag, "x": xx})
    res = serial_block(y, X, names, lagged_dep_idx=(1,))
    res["median"] = slim(
        serial_block(y, X, names, lagged_dep_idx=(1,), ar_order=2, wallis="median")
    )
    return {
        "formula": "y ~ ylag1 + x",
        "T": T,
        "seed": 20260106,
        "lagged_dep_vars": ["ylag1"],
        "results": res,
    }


@scenario
def s7_hetero_log2(out: Path) -> dict[str, Any]:
    """GP24 eq. (9) with v(x) = [log x]^2, z = (x, [log x]^2, x^2, cos x, cos 2x), T = 200."""
    seed = 20260107
    sim = hetero_dgp(200, skedastic="log2", random_state=seed)
    x = sim.X[:, 1]
    Z, zn = gp24_z(x)
    xn = ["(Intercept)", "x"]
    write_csv(out / "s7_data.csv", {"y": sim.y, "x": x})
    res = hetero_block(sim.y, sim.X, Z, zn, seed)
    res["one_se"] = slim(hetero_block(sim.y, sim.X, Z, zn, seed, lambda_rule="1se"))
    res["wls_s1"] = slim(hetero_block(sim.y, sim.X, Z[:, [0, 2]], ["x", "x^2"], seed, method="wls_s1"))
    res["wls_s2"] = slim(hetero_block(sim.y, sim.X, Z, zn, seed, method="wls_s2"))
    res["svr"] = slim(hetero_block(sim.y, sim.X, Z, zn, seed, method="svr"))
    res["ols"] = slim(hetero_block(sim.y, sim.X, Z, zn, seed, method="none"))
    res["public"] = lasso_public(sim.y, sim.X, xn, Z, zn, seed)
    return {"formula": "y ~ x", "T": 200, "seed": seed, "z": "gp24", "results": res}


@scenario
def s8_hetero_many(out: Path) -> dict[str, Any]:
    """GP24 Case 3 style: v = exp(0.2x + 0.2x^2) with 60 candidate covariates, T = 200."""
    seed = 20260108
    sim = hetero_dgp(200, skedastic="exp", n_extra=19, extra_corr=0.5, random_state=seed)
    x = sim.X[:, 1]
    W = sim.info["W"]
    cols = [x, x**2, np.log(x) ** 2]
    names = ["x", "x^2", "log(x)^2"]
    for j in range(W.shape[1]):
        w = W[:, j]
        cols += [w, np.cos(w), np.cos(2.0 * w)]
        names += [f"w{j + 1}", f"cos(w{j + 1})", f"cos(2w{j + 1})"]
    Z = np.column_stack(cols)
    data = {"y": sim.y, "x": x}
    data.update({f"w{j + 1}": W[:, j] for j in range(W.shape[1])})
    write_csv(out / "s8_data.csv", data)
    res = hetero_block(sim.y, sim.X, Z, names, seed)
    return {"formula": "y ~ x", "T": 200, "seed": seed, "z": "many", "results": res}


@scenario
def s9_combined(out: Path) -> dict[str, Any]:
    """Supplement S-4: AR(2) (1.34, -0.42) errors with v(x) = exp(0.2x + 0.2x^2), T = 200."""
    seed = 20260109
    sim = combined_dgp(200, ar=(1.34, -0.42), skedastic="exp", driver="x", random_state=seed)
    xn = ["(Intercept)", "x"]
    write_csv(out / "s9_data.csv", {"y": sim.y, "x": sim.X[:, 1]})
    res = combined_block(sim.y, sim.X, xn, seed)
    nofs = combined_block(sim.y, sim.X, xn, seed, first_stage_correction=False)
    res["no_first_stage"] = slim(nofs) | {"vcov_FCH": nofs["vcov_FCH"], "vcov_FC": nofs["vcov_FC"]}
    return {"formula": "y ~ x", "T": 200, "seed": seed, "z": "auto", "results": res}


@scenario
def s10_homoskedastic(out: Path) -> dict[str, Any]:
    """Homoskedastic white-noise errors: no loss (Lasso selects nothing, F-CH = F-C = OLS)."""
    from pyfgls.hetero.features import make_dictionary

    seed = 20260110
    sim = hetero_dgp(200, skedastic="homoskedastic", random_state=seed)
    xn = ["(Intercept)", "x"]
    Z, zn = make_dictionary(sim.X, xn)
    write_csv(out / "s10_data.csv", {"y": sim.y, "x": sim.X[:, 1]})
    res = hetero_block(sim.y, sim.X, Z, list(zn), seed)
    res["dictionary_names"] = list(zn)
    res["combined"] = combined_block(sim.y, sim.X, xn, seed)
    return {"formula": "y ~ x", "T": 200, "seed": seed, "z": "auto", "results": res}


# ------------------------------------------------------------------------- drift
#: Outputs of iterative numerical solvers (the MA(1) maximum likelihood behind GMA,
#: libsvm's SVR) depend on the BLAS kernels in the last bits, which the solver
#: amplifies. Across OpenBLAS's kernels (OPENBLAS_CORETYPE, i.e. across CPUs) they move
#: by up to ~3e-8 (gma) and ~1e-6 (svr cv_mse); every other output by less than 1e-10.
#: The R tests compare them at 1e-4 (gma_theta) or not at all (SVR CV errors).
SOLVER_BLOCKS = frozenset({"gma", "svr"})
SOLVER_RTOL = 1e-5


def _numbers(x: Any) -> list[float] | None:
    """The numbers of a (nested) numeric list or scalar; None if it holds anything else."""
    if isinstance(x, list):
        out: list[float] = []
        for v in x:
            f = _numbers(v)
            if f is None:
                return None
            out.extend(f)
        return out
    if isinstance(x, int | float) and not isinstance(x, bool):
        return [float(x)]
    return None


def _compare_numbers(
    a: list[float], b: list[float], path: str, rtol: float, out: list[str]
) -> None:
    """Norm-wise relative difference max|a - b| / max(|a|, |b|), as in the R tests."""
    if not a:
        return
    diff = [abs(x - y) for x, y in zip(a, b, strict=True)]
    i = max(range(len(diff)), key=diff.__getitem__)
    scale = max(max(map(abs, a)), max(map(abs, b)), 1e-300)
    if diff[i] > rtol * scale and diff[i] > 1e-14:
        where = f"[{i}]" if len(a) > 1 else ""
        out.append(f"{path}{where}: {a[i]!r} != {b[i]!r} (norm-wise {diff[i] / scale:.1e})")


def _compare(a: Any, b: Any, path: str, rtol: float, out: list[str]) -> None:
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f"{path}.{k}: present in only one version")
            else:
                sub = max(rtol, SOLVER_RTOL) if k in SOLVER_BLOCKS else rtol
                _compare(a[k], b[k], f"{path}.{k}", sub, out)
        return
    if isinstance(a, list) and isinstance(b, list) and len(a) != len(b):
        out.append(f"{path}: length {len(a)} != {len(b)}")
        return
    na, nb = _numbers(a), _numbers(b)
    if na is not None and nb is not None and len(na) == len(nb):
        _compare_numbers(na, nb, path, rtol, out)
    elif isinstance(a, list) and isinstance(b, list):
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            _compare(x, y, f"{path}[{i}]", rtol, out)
    elif a != b:
        out.append(f"{path}: {a!r} != {b!r}")


def check_drift(committed: Path, fresh: Path, rtol: float) -> int:
    """Compare regenerated fixtures with the committed ones; non-zero exit on drift."""
    problems: list[str] = []
    for f in sorted(fresh.glob("*.json")):
        if f.name == "MANIFEST.json":
            continue
        old = committed / f.name
        if not old.exists():
            problems.append(f"{f.name}: not committed")
            continue
        _compare(json.loads(old.read_text()), json.loads(f.read_text()), f.stem, rtol, problems)
    for f in sorted(fresh.glob("*.csv")):
        old = committed / f.name
        if not old.exists() or old.read_bytes() != f.read_bytes():
            problems.append(f"{f.name}: data differ")
    if problems:
        print(f"Golden fixtures drifted ({len(problems)} differences; first 50):")
        print("\n".join(problems[:50]))
        return 1
    print("No drift: regenerated fixtures agree with the committed ones.")
    return 0


# ------------------------------------------------------------------------------ main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    ap.add_argument("--out", default=str(root / "tests" / "testthat" / "fixtures" / "golden"))
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument(
        "--check",
        action="store_true",
        help="regenerate into a temporary directory and compare with the committed fixtures",
    )
    ap.add_argument(
        "--rtol", type=float, default=1e-10, help="norm-wise relative tolerance per array for --check"
    )
    args = ap.parse_args(argv)
    if args.check:
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            main(["--out", tmp] + (["--only", *args.only] if args.only else []))
            return check_drift(Path(args.out), Path(tmp), args.rtol)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if pyfgls.__version__ != PYFGLS_PIN:
        print(f"warning: pyfgls {pyfgls.__version__} != pinned {PYFGLS_PIN}", file=sys.stderr)

    manifest_path = out / "MANIFEST.json"
    manifest: dict[str, Any] = {}
    if manifest_path.exists() and args.only:
        manifest = json.loads(manifest_path.read_text())
    import scipy
    import sklearn

    manifest.update(
        {
            "pyfgls": pyfgls.__version__,
            "pyfgls_pin": PYFGLS_PIN,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "scikit-learn": sklearn.__version__,
            "python": platform.python_version(),
        }
    )
    files = manifest.setdefault("scenarios", {})
    for name, fn in SCENARIOS.items():
        if args.only and name not in args.only:
            continue
        spec = fn(out)
        spec["scenario"] = name
        spec["description"] = (fn.__doc__ or "").strip()
        path = out / f"{name.split('_')[0]}_golden.json"
        path.write_text(json.dumps(jsonable(spec), indent=None, separators=(",", ":")))
        data = out / f"{name.split('_')[0]}_data.csv"
        files[name] = {
            "json": path.name,
            "data": data.name,
            "sha256_data": hashlib.sha256(data.read_bytes()).hexdigest(),
        }
        r = spec["results"]
        print(f"{name}: k*={r.get('kstar')} selected={r.get('nonzero_set')} psi={r.get('psi_star')}")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    total = sum(p.stat().st_size for p in out.iterdir())
    print(f"total fixture size: {total / 1024:.1f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
