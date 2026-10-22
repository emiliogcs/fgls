#!/usr/bin/env python3
"""pyfgls side of the cross-language validation (see tools/compare_python_r.R).

Reads the data sets and cross-validation folds written by
``Rscript tools/compare_python_r.R simulate`` and fits them with pyfgls, using the
*same folds* (``foldid=`` when pyfgls supports it, otherwise the fold generator of
the pinned version is replaced for the duration of each fit).

    python tools/make_validation.py validation/cache [workers]
"""

from __future__ import annotations

import inspect
import json
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import importlib

import pyfgls

# the package re-exports the function adaptive_lasso under the submodule's name
al = importlib.import_module("pyfgls.hetero.adaptive_lasso")


def paper_z(x: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(
        np.column_stack([x, np.log(x) ** 2, x**2, np.cos(x), np.cos(2 * x)]),
        columns=["x", "log(x)^2", "x^2", "cos(x)", "cos(2x)"],
    )


def paper_w(w: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(
        np.column_stack([w, w**2, np.cos(w), np.cos(2 * w)]),
        columns=["w", "w^2", "cos(w)", "cos(2w)"],
    )


HAS_FOLDID = "foldid" in inspect.signature(pyfgls.FGLS.__init__).parameters


def with_folds(foldid: list[int] | None, make: Any) -> Any:
    """Build a model with the given folds (foldid= if available, else patched folds)."""
    if foldid is None:
        return make({})
    folds = np.asarray(foldid, dtype=np.intp) - 1
    if HAS_FOLDID:
        return make({"foldid": folds})
    original = al.make_folds
    al.make_folds = lambda n, K, scheme="random", random_state=None: folds.copy()
    try:
        return make({})
    finally:
        al.make_folds = original


def fit_one(kind: str, ds: dict[str, Any]) -> dict[str, Any]:
    y = pd.Series(np.asarray(ds["y"], dtype=float), name="y")
    x = np.asarray(ds["x"], dtype=float)
    X = pd.DataFrame({"(Intercept)": np.ones(x.size), "x": x})
    foldid = ds.get("foldid")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if kind in ("serial", "ma1"):
            res = pyfgls.SerialFGLS(y, X).fit()
            out: dict[str, Any] = {
                "beta": float(res.params.iloc[1]),
                "se": float(res.bse.iloc[1]),
                "kstar": int(res.ar_order),
                "rho1": float(res.rho.iloc[0]) if res.ar_order else None,
            }
            if kind == "serial":
                naive = pyfgls.SerialFGLS(y, X).fit(cov_type="naive")
                out["se_naive"] = float(naive.bse.iloc[1])
                out["share"] = float(res.nonexogeneity_share)
            else:
                g = pyfgls.SerialFGLS(y, X, ma1=True).fit()
                out["theta"] = float(g.serial.extra["MA coef theta"])
                out["beta_gma"] = float(g.params.iloc[1])
                out["se_gma"] = float(g.bse.iloc[1])
            return out
        if kind == "hetero":
            res = with_folds(
                foldid, lambda kw: pyfgls.LassoFGLS(y, X, z=paper_z(x), **kw).fit()
            )
            return {
                "beta": float(res.params.iloc[1]),
                "se": float(res.bse.iloc[1]),
                "psi": float(res.hetero.tuning["psi"]),
                "lambda": float(res.hetero.tuning["lambda"]),
                "df_hat": float(res.df_hat),
                "selected": "|".join(res.selected_z),
            }
        w = np.asarray(ds["w"], dtype=float)
        res = with_folds(foldid, lambda kw: pyfgls.FGLS(y, X, z=paper_w(w), **kw).fit())
        fc = res.stage_serial
        return {
            "beta": float(res.params.iloc[1]),
            "se": float(res.bse.iloc[1]),
            "kstar": int(res.ar_order),
            "beta_fc": float(fc.params.iloc[1]),
            "se_fc": float(fc.bse.iloc[1]),
            "psi": float(res.hetero.tuning["psi"]),
            "df_hat": float(res.df_hat),
            "selected": "|".join(res.selected_z),
        }


def _work(args: tuple[str, dict[str, Any]]) -> dict[str, Any]:
    return fit_one(*args)


def main(argv: list[str]) -> int:
    cache = Path(argv[0] if argv else "validation/cache")
    workers = int(argv[1]) if len(argv) > 1 else 1
    for f in sorted(cache.glob("*.json")):
        if f.name.endswith("_py.json"):
            continue
        spec = json.loads(f.read_text())
        jobs = [(spec["kind"], ds) for ds in spec["datasets"]]
        if workers > 1:
            with Pool(workers) as pool:
                rows = pool.map(_work, jobs, chunksize=8)
        else:
            rows = [_work(j) for j in jobs]
        keys = list(rows[0])
        out = {
            "design": spec["design"],
            "pyfgls": pyfgls.__version__,
            "foldid_mode": "foldid" if HAS_FOLDID else "patched make_folds",
            "results": {k: [r[k] for r in rows] for k in keys},
        }
        (cache / f"{spec['design']}_py.json").write_text(json.dumps(out))
        print(f"{spec['design']}: {len(rows)} fits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
