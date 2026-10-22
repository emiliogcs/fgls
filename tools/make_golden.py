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
from pyfgls.simulate import serial_dgp

PYFGLS_PIN = "0.1.0"
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


def public_api_check(y, X, names, **kw) -> dict[str, Any]:
    """End-to-end pyfgls result through the public API (sanity check of the internals)."""
    import pandas as pd

    df = pd.DataFrame(X, columns=names)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = pyfgls.SerialFGLS(pd.Series(y, name="y"), df, **kw).fit()
    return {"params": res.params.to_numpy(), "bse": res.bse.to_numpy(), "ar_order": res.ar_order}


# -------------------------------------------------------------------------- scenarios
@scenario
def s1_ar1_nonexog(out: Path) -> dict[str, Any]:
    """AR(1) errors (0.8), non-exogenous regressor (gamma = 0.5), PG26 Table-1 design."""
    sim = serial_dgp(500, ar=(0.8,), rho_x=0.0, gamma=0.5, random_state=20260101)
    names = ["(Intercept)", "x"]
    write_csv(out / "s1_data.csv", {"y": sim.y, "x": sim.X[:, 1]})
    res = serial_block(sim.y, sim.X, names)
    res["fixed_order_3"] = serial_block(sim.y, sim.X, names, ar_order=3)
    res["public"] = public_api_check(sim.y, sim.X, names)
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
    res["median"] = serial_block(y, X, names, lagged_dep_idx=(1,), ar_order=2, wallis="median")
    return {
        "formula": "y ~ ylag1 + x",
        "T": T,
        "seed": 20260106,
        "lagged_dep_vars": ["ylag1"],
        "results": res,
    }


# ------------------------------------------------------------------------------ main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    ap.add_argument("--out", default=str(root / "tests" / "testthat" / "fixtures" / "golden"))
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args(argv)
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
        print(f"{name}: k*={spec['results'].get('kstar')}")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    total = sum(p.stat().st_size for p in out.iterdir())
    print(f"total fixture size: {total / 1024:.1f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
