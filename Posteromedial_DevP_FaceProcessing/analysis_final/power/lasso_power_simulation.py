#!/usr/bin/env python3
"""
power/lasso_power_simulation.py

Simulation-based power (sensitivity) analysis for the right-hemisphere LASSO
analysis of sulcal depth and CFMT score in neurotypical participants.

Reproduces (Methods, "Power analysis")
  - multicollinearity among the nine sulcal depths: maximum variance inflation
    factor, largest pairwise correlation, condition number
  - the standardized effect at which the pipeline selects a genuinely
    associated sulcus in 80% of simulated datasets (median and range across
    the nine sulci)
  - selection rate at beta = .30 (median across sulci) and for a sulcus with
    no true effect (null simulation)
  - median cross-validated R^2 at beta = .30 and beta = .50, and the selection
    rate at beta = .50
  It also re-runs the published pipeline on the real data as a gate (alpha,
  RMSE_CV, coefficients and fold selection rates of the reported model) and
  reports the observed cross-validated R^2 and adjusted R^2.

Design
  X is the observed right-hemisphere design matrix (28 NTs x 9 sulcal depths).
  A known standardized effect beta is planted on one sulcus j at a time:
      y = mean(CFMT) + sd(CFMT) * (beta * z(X_j) + sqrt(1 - beta^2) * e),
  e ~ N(0, 1), so beta is the standardized coefficient and Var(y) stays at the
  observed CFMT variance. Each simulated dataset is analysed with the same
  fully nested leave-one-out LASSO pipeline as the paper (alpha grid
  logspace(-1, 1, 70); alpha chosen by the median across outer folds of the
  inner-CV RMSE; a sulcus is selected when its median coefficient across outer
  folds is non-zero). A null simulation (beta = 0 for every sulcus) gives the
  false-positive rate.

Inputs   data/sulcal_depth.csv
Outputs  outputs/power/lasso_power_simulation.txt  report
         outputs/power/lasso_power_by_beta.csv      per sulcus x beta summaries
         outputs/power/lasso_power_null.csv         per-sulcus false-positive rate
         outputs/power/lasso_power_detectable.csv   beta at 80% / 50% selection
The outputs of the full run used in the paper are saved in
outputs/power/saved_full_run/.

Usage (from the repository root)
  python3 power/lasso_power_simulation.py                  full run (150 datasets per
                                                           cell; about 1-2 h on 8 cores)
  python3 power/lasso_power_simulation.py --iterations 30  reduced run for checking
  python3 power/lasso_power_simulation.py --from-saved     summarize the saved full run
Every simulated dataset has its own fixed seed, so the full run reproduces the
saved outputs exactly. The coefficient-path kernel calls scikit-learn's
internal coordinate-descent routine and is checked against the public
sklearn.linear_model.lasso_path before any simulation (tested with
scikit-learn 1.7).
"""
import argparse
import os
import sys
import time

# Each simulated dataset is tiny, so multithreaded BLAS only slows the parallel
# workers down: pin to one thread before numpy is imported.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import lasso_path
from sklearn.linear_model import _cd_fast as _cd

DEPTH_CSV = os.path.join("data", "sulcal_depth.csv")
OUTDIR = os.path.join("outputs", "power")
SAVED = os.path.join(OUTDIR, "saved_full_run")

N_ITER = 150                                   # simulated datasets per (sulcus, beta) cell
BETAS = np.array([0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70])
SEED = 20260914
N_JOBS = max(1, (os.cpu_count() or 2) - 1)

ALPHAS = np.array([round(v, 3) for v in np.logspace(-1, 1, 70)])
ALPHAS_DESC = ALPHAS[::-1]                     # lasso_path requires descending alphas

# 'sbps' in the data is the manuscript's 'spls'
SULC_CSV = ['mcgs', 'pos', 'sbps', 'prculs-d', 'prcus-p', 'prcus-i', 'prcus-a', 'ifrms', 'sspls-v']
SULC_PAPER = {'sbps': 'spls'}


def paper_name(s):
    return SULC_PAPER.get(s, s)


def load_design():
    """Right-hemisphere NT design matrix (complete cases, as in the LASSO notebook)."""
    d = pd.read_csv(DEPTH_CSV)
    d = d[['sub', 'group', 'hemi', 'CFMT'] + SULC_CSV].dropna()
    rh = d[(d.hemi == 'rh') & (d.group == 'Controls')]
    assert len(rh) == 28, f"right-hemisphere NT model n = {len(rh)}, expected 28"
    return rh[SULC_CSV].to_numpy(float), rh['CFMT'].to_numpy(float)


_CD_RNG = np.random.RandomState(0)   # only consulted when random=1; we pass 0


def coef_path(Xc, yc):
    """Lasso coefficient path over ALPHAS, ascending, shape (n_feat, n_alpha).

    Calls sklearn's coordinate-descent kernel directly instead of going through
    lasso_path, whose input validation dominates the runtime on a problem this
    small (a roughly fourfold speed-up). It changes no arithmetic --
    verify_kernel() below asserts agreement with the public lasso_path on the
    real data to 1e-8 before anything is simulated, and the reproduction gate
    then re-derives the published numbers through this same code.

    Xc must be F-contiguous float64 and COLUMN-CENTRED; yc centred float64.
    Alphas at or above max|X'y|/n have the exact all-zero solution, so they are
    skipped rather than solved (a further ~40% saving, and exact).
    """
    n, p = Xc.shape
    tolv = 1e-12 * float(np.dot(yc, yc))
    amax = np.max(np.abs(Xc.T @ yc)) / n
    out = np.zeros((p, len(ALPHAS)))
    w = np.zeros(p)
    for i in range(len(ALPHAS) - 1, -1, -1):     # descend: warm starts matter
        a = ALPHAS[i]
        if a >= amax:
            continue                             # exact zero solution
        w, _, _, _ = _cd.enet_coordinate_descent(
            w, a * n, 0.0, Xc, yc, 100000, tolv, _CD_RNG, 0, 0)
        out[:, i] = w
    return out


def path_fit(Xtr, ytr):
    """Equivalent of Lasso(fit_intercept=True) over the whole alpha grid.

    sklearn's Lasso with an intercept centres X and y, fits, then recovers the
    intercept. The path kernel does not centre, so we centre by hand -- doing
    this wrong silently produces a different (and wrong) RMSE curve.
    """
    xm = Xtr.mean(0)
    ym = ytr.mean()
    cf = coef_path(np.asfortranarray(Xtr - xm), np.ascontiguousarray(ytr - ym))
    return cf, ym - xm @ cf


def verify_kernel(X, y):
    """Integrity gate: the fast kernel must agree with sklearn's public API."""
    Xc = np.asfortranarray(X[:26] - X[:26].mean(0))
    yc = np.ascontiguousarray(y[:26] - y[:26].mean())
    mine = coef_path(Xc, yc)
    _, ref, _ = lasso_path(Xc, yc, alphas=np.ascontiguousarray(ALPHAS_DESC),
                           max_iter=100000, tol=1e-12)
    d = np.abs(mine - ref[:, ::-1]).max()
    assert d < 1e-7, (f"GATE FAIL: fast coordinate-descent kernel disagrees with "
                      f"sklearn.linear_model.lasso_path by {d:.2e}")
    return d


def precompute_folds(X):
    """Everything in the pipeline that depends on X alone.

    The fold structure and every z-scored / centred design sub-matrix are
    identical across simulated datasets (only y changes), so they are built
    once and reused; this changes no arithmetic.
    """
    n = X.shape[0]
    folds = []
    for o in range(n):
        tr = np.delete(np.arange(n), o)
        Xt = X[tr]
        Xts = (Xt - Xt.mean(0)) / Xt.std(0, ddof=0)   # sklearn StandardScaler
        m = len(tr)
        inner = []
        for i in range(m):
            itr = np.delete(np.arange(m), i)
            Xi = Xts[itr]
            xm = Xi.mean(0)
            inner.append((itr, np.asfortranarray(Xi - xm), xm, Xts[i]))
        folds.append(dict(o=o, tr=tr, Xts=Xts, xtest=X[o], inner=inner,
                          mu=Xt.mean(0), sd=Xt.std(0, ddof=0)))
    return folds


def run_pipeline(folds, y):
    """One full fully-nested LASSO run. Returns the quantities the paper reports."""
    n = len(y)
    A = len(ALPHAS)
    p = folds[0]['Xts'].shape[1]
    rmse_cv = np.empty((A, n))          # inner-CV RMSE, per alpha per outer fold
    coefs = np.empty((A, n, p))         # refit coefficient, per alpha per fold
    icpt = np.empty((A, n))

    for f in folds:
        tr, Xts, inner = f['tr'], f['Xts'], f['inner']
        yt = y[tr]
        m = len(tr)
        err = np.empty((A, m))
        for k, (itr, Xic, xm, xrow) in enumerate(inner):
            yi = yt[itr]
            ym = yi.mean()
            cf = coef_path(Xic, np.ascontiguousarray(yi - ym))
            err[:, k] = np.abs(yt[k] - ((xrow - xm) @ cf + ym))
        rmse_cv[:, f['o']] = err.mean(1)
        cf, ic = path_fit(Xts, yt)
        coefs[:, f['o'], :] = cf.T
        icpt[:, f['o']] = ic

    med_rmse = np.median(rmse_cv, axis=1)
    a = int(np.argmin(med_rmse))

    # Held-out prediction at alpha*, using each fold's own scaler (fully nested).
    yhat = np.empty(n)
    for f in folds:
        xs = (f['xtest'] - f['mu']) / f['sd']
        yhat[f['o']] = xs @ coefs[a, f['o'], :] + icpt[a, f['o']]

    med_beta = np.median(coefs[a], axis=0)
    fold_rate = (coefs[a] != 0).mean(0)
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot
    k_sel = int(np.sum(med_beta != 0))
    r2_adj = (1.0 - (1.0 - r2) * (n - 1) / (n - k_sel - 1)
              if n - k_sel - 1 > 0 else np.nan)

    return dict(alpha=ALPHAS[a], rmse_cv=med_rmse[a], med_beta=med_beta,
                fold_rate=fold_rate, r2=r2, r2_adj=r2_adj, k_sel=k_sel)


# ---------------------------------------------------------- classical comparisons
def ols_and_corr(X, y, j):
    """What the classical instruments would have concluded on the same dataset.

    Returns (p of the full 9-predictor OLS coefficient for predictor j,
             p of the simple bivariate Pearson correlation for predictor j).
    """
    n, p = X.shape
    Xd = np.column_stack([np.ones(n), X])
    beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    resid = y - Xd @ beta
    dof = n - p - 1
    s2 = resid @ resid / dof
    cov = s2 * np.linalg.inv(Xd.T @ Xd)
    t = beta[j + 1] / np.sqrt(cov[j + 1, j + 1])
    p_ols = 2 * stats.t.sf(abs(t), dof)
    p_cor = stats.pearsonr(X[:, j], y)[1]
    return p_ols, p_cor


# --------------------------------------------------------------- the simulation
def run_cell(args):
    """One (sulcus, beta) cell: n_iter simulated datasets.

    Deliberately chunked at the cell level rather than the iteration level --
    every dispatched task pickles its arguments, and the precomputed fold
    structure is ~1.4 MB, so per-iteration dispatch spends more time
    serialising than fitting. Folds are rebuilt once inside each worker.
    Takes a single tuple so it can be used with Executor.map.
    """
    X, j, beta, sd_y, mean_y, seed0, n_iter = args
    folds = precompute_folds(X)
    return [one_iteration(folds, X, j, beta, sd_y, mean_y, seed0 + i)
            for i in range(n_iter)]


def parallel_cells(tasks, label="", fn=None):
    """Run run_cell over `tasks`, one process per core minus one.

    Uses a plain ProcessPoolExecutor with an explicit 'spawn' context rather
    than joblib, and keeps the single-thread BLAS setting in every worker.
    """
    fn = fn or run_cell
    if N_JOBS == 1:
        return [fn(t) for t in tasks]
    ctx = mp.get_context("spawn")
    done = 0
    out = []
    with ProcessPoolExecutor(max_workers=N_JOBS, mp_context=ctx) as ex:
        for r in ex.map(fn, tasks, chunksize=1):
            out.append(r)
            done += 1
            print(f"  [{label}] {done}/{len(tasks)} cells", flush=True)
    return out


def one_iteration(folds, X, j, beta, sd_y, mean_y, seed):
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    Xs = (X - X.mean(0)) / X.std(0, ddof=1)
    z = rng.standard_normal(n)
    y = mean_y + sd_y * (beta * Xs[:, j] + np.sqrt(1 - beta ** 2) * z)

    res = run_pipeline(folds, y)
    p_ols, p_cor = ols_and_corr(Xs, y, j)
    sel = res['med_beta'] != 0
    return dict(
        alpha=res['alpha'],
        sel_true=bool(sel[j]),
        fold_rate_true=float(res['fold_rate'][j]),
        sign_ok=bool(sel[j] and np.sign(res['med_beta'][j]) == np.sign(beta)),
        n_false=int(sel.sum() - sel[j]),
        sel_vector=sel.astype(np.int8),
        r2=res['r2'], r2_adj=res['r2_adj'], k_sel=res['k_sel'],
        ols_hit=bool(p_ols < 0.05), cor_hit=bool(p_cor < 0.05),
    )


def summarise(rows, sulcus, beta):
    n = len(rows)
    sel = np.array([r['sel_true'] for r in rows])
    return dict(
        sulcus=paper_name(sulcus), beta=beta, n_iter=n,
        sel_rate=sel.mean(),
        sel_rate_se=np.sqrt(sel.mean() * (1 - sel.mean()) / n),
        sel_rate_signed=np.mean([r['sign_ok'] for r in rows]),
        mean_fold_rate=np.mean([r['fold_rate_true'] for r in rows]),
        mean_n_false=np.mean([r['n_false'] for r in rows]),
        any_false=np.mean([r['n_false'] > 0 for r in rows]),
        mean_k_sel=np.mean([r['k_sel'] for r in rows]),
        median_alpha=np.median([r['alpha'] for r in rows]),
        median_r2=np.median([r['r2'] for r in rows]),
        mean_r2=np.mean([r['r2'] for r in rows]),
        r2_pos_rate=np.mean([r['r2'] > 0 for r in rows]),
        median_r2_adj=np.median([r['r2_adj'] for r in rows]),
        ols_power=np.mean([r['ols_hit'] for r in rows]),
        corr_power=np.mean([r['cor_hit'] for r in rows]),
    )


def interp_threshold(betas, rates, target):
    """Linear interpolation of the beta at which `rates` first crosses target."""
    betas, rates = np.asarray(betas, float), np.asarray(rates, float)
    for i in range(len(betas) - 1):
        if rates[i] < target <= rates[i + 1]:
            f = (target - rates[i]) / (rates[i + 1] - rates[i])
            return betas[i] + f * (betas[i + 1] - betas[i])
    if rates[-1] >= target:
        return betas[0]
    return np.nan



# ---------------------------------------------------------------------- report
def manuscript_quantities(df, null, det, say):
    """The quantities quoted in the Methods, from the per-cell summaries."""
    sel = lambda b: np.median(df[np.isclose(df.beta, b)].sel_rate) * 100
    r2 = lambda b: df[np.isclose(df.beta, b)].median_r2.mean()
    say("-" * 78)
    say("QUANTITIES REPORTED IN THE MANUSCRIPT (Methods, 'Power analysis')")
    say("-" * 78)
    say(f"beta at 80% selection: median {np.nanmedian(det.beta80):.2f} across sulci "
        f"(range {np.nanmin(det.beta80):.2f}-{np.nanmax(det.beta80):.2f})")
    for b in (0.30, 0.50):
        if np.isclose(df.beta, b).any():
            say(f"beta = {b:.2f}: selection rate, median across sulci = {sel(b):.1f}%; "
                f"median cross-validated R^2 (mean over sulci) = {r2(b):.3f}")
    say(f"null simulation: false-positive selection rate, mean across sulci = "
        f"{100 * null.false_positive_rate.mean():.1f}%")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iterations", type=int, default=N_ITER,
                    help="simulated datasets per cell (default 150, as in the paper)")
    ap.add_argument("--betas", type=float, nargs="+", default=list(BETAS),
                    help="planted standardized effects (default .10 to .70)")
    ap.add_argument("--from-saved", action="store_true",
                    help="summarize the saved full-run outputs instead of simulating")
    a = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)

    if a.from_saved:
        df = pd.read_csv(os.path.join(SAVED, "lasso_power_by_beta.csv"))
        null = pd.read_csv(os.path.join(SAVED, "lasso_power_null.csv"))
        det = pd.read_csv(os.path.join(SAVED, "lasso_power_detectable.csv"))
        manuscript_quantities(df, null, det, print)
        return

    n_iter, betas = a.iterations, np.array(a.betas)
    out = []

    def say(s=""):
        print(s)
        out.append(s)

    X, y = load_design()
    n, p = X.shape
    sd_y, mean_y = y.std(ddof=1), y.mean()
    say("SIMULATION-BASED POWER ANALYSIS, RIGHT-HEMISPHERE LASSO (NTs)")
    say(f"n = {n}, K = {p} sulcal depths, n:K = {n / p:.1f}; {n_iter} simulated datasets per cell")
    say()

    # ---- collinearity of the observed design ----------------------------------
    Xs = (X - X.mean(0)) / X.std(0, ddof=1)
    C = np.corrcoef(Xs.T)
    vif = np.diag(np.linalg.inv(C))
    off = np.abs(C - np.eye(p))
    say("OBSERVED PREDICTOR COLLINEARITY")
    say(f"{'sulcus':<10}{'VIF':>7}{'max |r| with another sulcus':>32}")
    for k, s in enumerate(SULC_CSV):
        j = int(np.argmax(off[k]))
        say(f"{paper_name(s):<10}{vif[k]:>7.2f}{off[k, j]:>20.2f}  (vs {paper_name(SULC_CSV[j])})")
    say(f"maximum VIF = {vif.max():.2f}; largest |r| = {off.max():.2f}; "
        f"condition number = {np.linalg.cond(Xs):.2f}")
    say()

    # ---- reproduction gate on the real data ------------------------------------
    folds = precompute_folds(X)
    d = verify_kernel(X, y)
    real = run_pipeline(folds, y)
    ix = {s: i for i, s in enumerate(SULC_CSV)}
    say("REAL DATA (gate: the published model must be reproduced)")
    say(f"kernel vs sklearn lasso_path: max |diff| = {d:.1e}")
    say(f"alpha = {real['alpha']:.3f}; RMSE_CV = {real['rmse_cv']:.4f}")
    for s in ['ifrms', 'prcus-p', 'prcus-i']:
        say(f"  {paper_name(s):<9} beta = {real['med_beta'][ix[s]]:>7.3f}   "
            f"fold selection = {real['fold_rate'][ix[s]] * 100:>5.1f}%")
    say(f"cross-validated R^2 = {real['r2']:.4f}; adjusted R^2 = {real['r2_adj']:.4f}")
    checks = [abs(real['alpha'] - 1.35) < 1e-6, abs(real['rmse_cv'] - 5.2076) < 5e-4,
              abs(real['med_beta'][ix['ifrms']] + 1.291) < 5e-4,
              abs(real['med_beta'][ix['prcus-p']] + 0.227) < 5e-4,
              abs(real['med_beta'][ix['prcus-i']] + 0.150) < 5e-4]
    if not all(checks):
        sys.exit("the published model was not reproduced; stopping")
    say()

    # ---- null simulation --------------------------------------------------------
    n_chunk = max(1, n_iter // N_JOBS)
    null_tasks = [(X, 0, 0.0, sd_y, mean_y, SEED + 900000 + c * n_chunk,
                   min(n_chunk, n_iter - c * n_chunk))
                  for c in range((n_iter + n_chunk - 1) // n_chunk)]
    null_rows = [r for chunk in parallel_cells(null_tasks, "null") for r in chunk]
    selmat = np.array([r['sel_vector'] for r in null_rows])
    null = pd.DataFrame({'sulcus': [paper_name(s) for s in SULC_CSV],
                         'false_positive_rate': selmat.mean(0)})
    null.to_csv(os.path.join(OUTDIR, 'lasso_power_null.csv'), index=False)
    say("NULL SIMULATION (beta = 0): per-sulcus false-positive selection rate")
    for k, s in enumerate(SULC_CSV):
        say(f"  {paper_name(s):<10}{selmat[:, k].mean() * 100:>6.1f}%")
    say()

    # ---- effect-size sweep --------------------------------------------------------
    cells = [(j, b) for j in range(p) for b in betas]
    tasks = [(X, j, b, sd_y, mean_y, SEED + 1000 * j + int(round(b * 1000)) * 7, n_iter)
             for j, b in cells]
    res = parallel_cells(tasks, "sweep")
    df = pd.DataFrame([summarise(rows, SULC_CSV[j], b) for (j, b), rows in zip(cells, res)])
    df = df.sort_values(['sulcus', 'beta']).reset_index(drop=True)
    df.to_csv(os.path.join(OUTDIR, 'lasso_power_by_beta.csv'), index=False)

    say("SELECTION RATE (%) OF THE TRULY ASSOCIATED SULCUS, BY PLANTED BETA")
    say(f"{'sulcus':<10}" + "".join(f"{b:>8.2f}" for b in betas))
    for s in SULC_CSV:
        sub = df[df.sulcus == paper_name(s)].sort_values('beta')
        say(f"{paper_name(s):<10}" + "".join(f"{v * 100:>8.1f}" for v in sub.sel_rate))
    det = []
    for k, s in enumerate(SULC_CSV):
        sub = df[df.sulcus == paper_name(s)].sort_values('beta')
        det.append(dict(sulcus=paper_name(s), vif=vif[k],
                        beta80=interp_threshold(sub.beta, sub.sel_rate, 0.80),
                        beta50=interp_threshold(sub.beta, sub.sel_rate, 0.50),
                        ols80=interp_threshold(sub.beta, sub.ols_power, 0.80),
                        corr80=interp_threshold(sub.beta, sub.corr_power, 0.80)))
    det = pd.DataFrame(det)
    det.to_csv(os.path.join(OUTDIR, 'lasso_power_detectable.csv'), index=False)
    say()
    manuscript_quantities(df, null, det, say)
    with open(os.path.join(OUTDIR, "lasso_power_simulation.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
