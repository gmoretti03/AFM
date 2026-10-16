"""Shared functions for notebooks 02-05 of the AFM pipeline.

Contents
--------
- reading of the segment0 (extend/approach) curves exported by notebook 01;
- robust linear baseline, noise, contact point estimate;
- Hertz-like fit  F(delta) = A * delta^(3/2)  with free contact point;
- fit variants (force window, threshold) for the systematic-uncertainty study;
- statistical utilities (bootstrap, quantiles) and a common plot style.

Conventions (as in notebook 01):
- ``z_tip_um`` DECREASES along the extend curve: the baseline (tip far from the
  surface) is at the START of the array, the contact towards the END.
- ``A`` is an *effective* elastic parameter in nN / um^(3/2): it is not the absolute
  Young modulus (tip radius and Poisson ratio are missing).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
CODE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = CODE_DIR / "afm_output"          # regenerable cache, not tracked
RESULTS_DIR = CODE_DIR / "results"            # light results, tracked by git
MANIFEST_PATH = OUTPUT_DIR / "manifests" / "segment0_manifest.csv"
FIT_TABLE_PATH = RESULTS_DIR / "fit_results_segment0.csv"

GROUP_ORDER = [
    "Petri20k_MCF7R_t48h_cito_20260225",
    "Petri20k_MCF7R_t48h_nucleo_20260225",
    "T47DR_t48h_cito_20260318",
    "T47DR_t48h_nucleo_20260318",
]


# ----------------------------------------------------------------------------
# Common plot style
# ----------------------------------------------------------------------------
def set_plot_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight",
        "axes.grid": True, "grid.alpha": 0.3, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False, "font.size": 10,
    })


# ----------------------------------------------------------------------------
# Analysis parameters
# ----------------------------------------------------------------------------
@dataclass(frozen=True)
class FitConfig:
    """Preprocessing and fit parameters. The variants feed the systematic study."""
    baseline_frac: float = 0.30      # initial fraction of the curve used for the baseline
    threshold_sigma: float = 8.0     # contact threshold: k * noise sigma
    threshold_fmax_frac: float = 0.05  # ... or fraction of the maximum force
    min_consecutive: int = 20        # consecutive points above threshold to declare contact
    f_low_frac: float = 0.08         # fit window: from f_low_frac * Fmax ...
    f_high_frac: float = 0.80        # ... a f_high_frac * Fmax
    min_fit_points: int = 80

    def tag(self) -> str:
        return (f"bl{self.baseline_frac:g}_thr{self.threshold_sigma:g}_"
                f"lo{self.f_low_frac:g}_hi{self.f_high_frac:g}")


NOMINAL = FitConfig()

# Variants for the systematic study (notebook 03): one change at a time.
VARIANTS = {
    "nominal": NOMINAL,
    "baseline_20%": FitConfig(baseline_frac=0.20),
    "baseline_40%": FitConfig(baseline_frac=0.40),
    "threshold_6sigma": FitConfig(threshold_sigma=6.0),
    "threshold_12sigma": FitConfig(threshold_sigma=12.0),
    "window_low_4%": FitConfig(f_low_frac=0.04),
    "window_low_15%": FitConfig(f_low_frac=0.15),
    "window_high_60%": FitConfig(f_high_frac=0.60),
    "window_high_95%": FitConfig(f_high_frac=0.95),
}


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------
def read_curve(csv_path) -> tuple[np.ndarray, np.ndarray]:
    """Return (z_um, force_nN) of a segment0 CSV, without NaN."""
    df = pd.read_csv(csv_path, usecols=["z_tip_um", "force_nN"])
    z = df["z_tip_um"].to_numpy(float)
    f = df["force_nN"].to_numpy(float)
    ok = np.isfinite(z) & np.isfinite(f)
    return z[ok], f[ok]


def load_manifest() -> pd.DataFrame:
    """Manifest of notebook 01 (needs afm_output/, regenerable with notebook 01)."""
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(
            f"{MANIFEST_PATH} does not exist: run notebook 01 on the raw data first."
        )
    return pd.read_csv(MANIFEST_PATH)


def load_fit_table() -> pd.DataFrame:
    """Fit table produced by notebook 02 (tracked in results/)."""
    return pd.read_csv(FIT_TABLE_PATH)


# ----------------------------------------------------------------------------
# Baseline, noise, contact
# ----------------------------------------------------------------------------
def robust_linear_baseline(z, f, frac, n_iter=5, clip=3.0):
    """Baseline line on the initial fraction of the curve, with iterative clipping.

    Returns (coefficients [slope, intercept], robust sigma).
    """
    n = max(50, int(frac * len(z)))
    zb, fb = z[:n], f[:n]
    keep = np.ones(n, bool)
    for _ in range(n_iter):
        coef = np.polyfit(zb[keep], fb[keep], 1)
        res = fb - np.polyval(coef, zb)
        sigma = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
        if not np.isfinite(sigma) or sigma == 0:
            sigma = np.std(res[keep])
        new_keep = np.abs(res) < clip * sigma
        if new_keep.sum() < 20 or np.array_equal(new_keep, keep):
            break
        keep = new_keep
    return coef, float(sigma)


def smooth(f, n):
    win = min(101, max(11, (n // 50) | 1))
    if win >= n:
        win = max(5, (n // 5) | 1)
    try:
        return savgol_filter(f, window_length=win, polyorder=2)
    except Exception:
        return f.copy()


def estimate_contact(z, fcorr, noise, cfg: FitConfig):
    """First point (from the start) where the smoothed force stays above threshold for
    ``min_consecutive`` consecutive points. Returns (index, f_smooth, threshold)."""
    n = len(z)
    f_s = smooth(fcorr, n)
    thr = max(cfg.threshold_sigma * noise, cfg.threshold_fmax_frac * np.nanmax(fcorr))
    above = f_s > thr
    k = cfg.min_consecutive
    run = np.convolve(above.astype(int), np.ones(k, int), mode="valid")
    hits = np.flatnonzero(run == k)
    if len(hits):
        return int(hits[0]), f_s, float(thr)
    return -1, f_s, float(thr)


# ----------------------------------------------------------------------------
# Hertz-like model and fit
# ----------------------------------------------------------------------------
def hertz_model(z_um, A, zc_um):
    """F[nN] = A * delta^(3/2), delta = zc - z (positive beyond the contact)."""
    delta = np.clip(zc_um - np.asarray(z_um), 0, None)
    return A * delta ** 1.5


def fit_curve(z, f, cfg: FitConfig = NOMINAL, return_arrays=False) -> dict:
    """Preprocessing + Hertz-like fit of an extend curve. Never raises:
    on problems it returns ``fit_status != 'ok'`` and the reason in ``error``."""
    out = {"n_points": int(len(z)), "fit_status": "failed", "error": ""}
    try:
        if len(z) < 200:
            raise ValueError("too few points")
        coef, noise = robust_linear_baseline(z, f, cfg.baseline_frac)
        fcorr = f - np.polyval(coef, z)
        fmax = float(np.nanmax(fcorr))
        out.update(baseline_slope_nN_per_um=float(coef[0]),
                   baseline_noise_pN=noise * 1e3, fmax_nN=fmax,
                   snr=fmax / noise if noise > 0 else np.nan)

        idx0, f_s, thr = estimate_contact(z, fcorr, noise, cfg)
        if idx0 < 0:
            raise ValueError("no contact above threshold")
        zc0 = float(z[idx0])
        out.update(contact_idx_initial=idx0, contact_z_initial_um=zc0, threshold_nN=thr,
                   jump_score=float(np.max(np.abs(np.diff(f_s))) / noise),
                   sat_fraction=float(np.mean(fcorr > 0.99 * fmax)))

        f_low = max(thr, cfg.f_low_frac * fmax)
        f_high = cfg.f_high_frac * fmax
        mask = (z < zc0) & (fcorr > f_low) & (fcorr < f_high)
        if mask.sum() < cfg.min_fit_points:
            mask = (z < zc0) & (fcorr > thr) & (fcorr < 0.90 * fmax)
        if mask.sum() < cfg.min_fit_points:
            raise ValueError("fit region too small")
        zf, ff = z[mask], fcorr[mask]

        d0 = max(zc0 - zf.min(), 1e-9)
        A0 = ff.max() / d0 ** 1.5
        zr = z.max() - z.min()
        popt, pcov = curve_fit(hertz_model, zf, ff, p0=[A0, zc0],
                               bounds=([0, zc0 - 0.35 * zr], [np.inf, zc0 + 0.20 * zr]),
                               maxfev=20000)
        A, zc = popt
        A_err, zc_err = (np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else (np.nan, np.nan))
        res = ff - hertz_model(zf, A, zc)
        rss, tss = float(np.sum(res ** 2)), float(np.sum((ff - ff.mean()) ** 2))
        r2 = 1 - rss / tss if tss > 0 else np.nan
        out.update(
            contact_z_fit_um=float(zc), contact_z_fit_err_um=float(zc_err),
            fit_force_low_nN=float(f_low), fit_force_high_nN=float(f_high),
            n_fit_points=int(mask.sum()),
            A_nN_per_um_3over2=float(A), A_err_nN_per_um_3over2=float(A_err),
            max_indentation_um=float(zc - zf.min()),
            contact_shift_rel=float((zc0 - zc) / zr),
            exponent_hat=_free_exponent(zf, ff, zc),
            fit_rmse_nN=float(np.sqrt(np.mean(res ** 2))), fit_rmse_rel=float(np.sqrt(np.mean(res ** 2)) / fmax),
            fit_r2=float(r2),
            fit_status="ok",
        )
        out["fit_quality_flag"] = quality_flag(out)
        if return_arrays:
            out["_arrays"] = dict(z=z, fcorr=fcorr, f_smooth=f_s, zf=zf, ff=ff,
                                  baseline=np.polyval(coef, z))
    except Exception as e:  # noqa: BLE001 - the batch fit must not stop
        out["error"] = str(e)
        out["fit_quality_flag"] = "failed"
    return out


def _free_exponent(zf, ff, zc) -> float:
    """Exponent p of F ~ delta^p from a log-log regression in the fit window
    (contact point fixed at the fitted value). Hertz sphere: 1.5; pyramidal tip: 2."""
    d = zc - zf
    m = (d > 0) & (ff > 0)
    if m.sum() < 20:
        return float("nan")
    return float(np.polyfit(np.log(d[m]), np.log(ff[m]), 1)[0])


def quality_flag(r: dict) -> str:
    """Fit quality flag (criteria of prototype 02_v4/03): good / suspicious / bad."""
    A, r2, n = r.get("A_nN_per_um_3over2", np.nan), r.get("fit_r2", np.nan), r.get("n_fit_points", 0)
    if not (np.isfinite(A) and A > 0 and np.isfinite(r2)):
        return "bad"
    if r2 > 0.90 and n >= 100:
        return "good"
    if r2 > 0.75 and n >= 80:
        return "suspicious"
    return "bad"


def fit_manifest(manifest: pd.DataFrame, cfg: FitConfig = NOMINAL, root: Path = CODE_DIR) -> pd.DataFrame:
    """Fit all curves of the manifest. Metadata columns are merged into the result."""
    rows = []
    for _, m in manifest.iterrows():
        z, f = read_curve(root / m["csv_path"])
        r = fit_curve(z, f, cfg)
        r.update(file_name=m["file_name"], group=m["group"], cell_line=m["cell_line"],
                 region=m["region"], date=m["date"])
        rows.append(r)
    df = pd.DataFrame(rows)
    df["logA"] = np.where(df["A_nN_per_um_3over2"] > 0, np.log(df["A_nN_per_um_3over2"]), np.nan)
    return df


# ----------------------------------------------------------------------------
# Basic statistics
# ----------------------------------------------------------------------------
def good_sample(fit: pd.DataFrame, flags=("good",)) -> pd.DataFrame:
    """Curves used in the statistical analysis: fit ok, A > 0, requested flags."""
    m = (fit["fit_status"] == "ok") & fit["fit_quality_flag"].isin(flags) & np.isfinite(fit["logA"])
    return fit[m].copy()


def bootstrap_stat(x, stat=np.median, n_boot=5000, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    return stat(x[idx], axis=1)


def short_group(g: str) -> str:
    """'Petri20k_MCF7R_t48h_cito_20260225' -> 'MCF7R cito' (cito = cytoplasm, nucleo = nucleus)."""
    toks = g.split("_")
    line = next(t for t in toks if t in ("MCF7R", "T47DR"))
    region = next(t for t in toks if t in ("cito", "nucleo"))
    return f"{line} {region}"


# ----------------------------------------------------------------------------
# Retract (segment1) features: adhesion and hysteresis
# ----------------------------------------------------------------------------
def retract_features(csv_path_seg0, root: Path = CODE_DIR, baseline_frac=0.30) -> dict:
    """Adhesion and work of adhesion from segment 1 (retract), if present.

    The retract baseline is at the END of the array (tip far from the surface).
    """
    p = Path(root) / str(csv_path_seg0).replace("_segment0", "_segment1")
    out = {"adhesion_nN": np.nan, "adhesion_work_nN_um": np.nan}
    if not p.exists():
        return out
    z, f = read_curve(p)
    if len(z) < 200:
        return out
    n = max(50, int(baseline_frac * len(z)))
    coef = np.polyfit(z[-n:], f[-n:], 1)
    fc = smooth(f - np.polyval(coef, z), len(z))
    out["adhesion_nN"] = float(max(0.0, -np.nanmin(fc)))
    neg = np.clip(-fc, 0, None)
    out["adhesion_work_nN_um"] = float(np.trapezoid(neg, np.abs(z)))
    return out


# ----------------------------------------------------------------------------
# Diagnostic plot of one curve
# ----------------------------------------------------------------------------
def plot_fit(ax, res: dict, title: str = ""):
    a = res["_arrays"]
    ax.plot(a["z"], a["fcorr"], lw=0.6, alpha=0.6, label="force - baseline")
    ax.plot(a["z"], a["f_smooth"], lw=1.0, label="smoothed")
    ax.scatter(a["zf"], a["ff"], s=3, alpha=0.4, color="C2", label="fit points")
    zg = np.linspace(a["zf"].min(), res["contact_z_fit_um"], 200)
    ax.plot(zg, hertz_model(zg, res["A_nN_per_um_3over2"], res["contact_z_fit_um"]),
            color="C3", lw=2, label=r"fit $A\,\delta^{3/2}$")
    ax.axvline(res["contact_z_fit_um"], ls="-.", color="k", lw=0.8, label="contact (fit)")
    ax.axhline(res["threshold_nN"], ls=":", color="gray", lw=0.8, label="threshold")
    ax.set_xlabel("vertical tip position [µm]")
    ax.set_ylabel("force [nN]")
    ax.set_title(title, fontsize=9)
