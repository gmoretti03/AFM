"""Regenerate the figures of the report from Code_General/results (and, for Fig. 1, from afm_output).

Run:  python make_figures.py      (kernel/env: environment.yml)
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

CODE = Path(__file__).resolve().parents[2] / "Code_General"
sys.path.insert(0, str(CODE))
import afm_tools as T

OUT = Path(__file__).resolve().parent / "figures"; OUT.mkdir(exist_ok=True)
T.set_plot_style(); plt.rcParams.update({"font.size": 8, "figure.dpi": 150})
R = CODE / "results"
fit = T.load_fit_table(); ana = T.good_sample(fit)
s3, s4, s5 = (json.load(open(R / f"0{i}_summary.json")) for i in (3, 4, 5))
lines = ["MCF7R", "T47DR"]

# Fig 1: example fit (needs the raw-data cache; skipped if absent)
if T.MANIFEST_PATH.exists():
    m = T.load_manifest(); row = m[m.group.str.contains("MCF7R_t48h_cito")].iloc[3]
    z, f = T.read_curve(CODE / row.csv_path); res = T.fit_curve(z, f, return_arrays=True)
    fig, ax = plt.subplots(figsize=(3.4, 2.4)); T.plot_fit(ax, res, ""); ax.legend(fontsize=5.5, loc="upper right")
    fig.savefig(OUT / "fit_example.pdf"); plt.close(fig)

# Fig 2: log A per group
fig, axs = plt.subplots(2, 2, figsize=(3.4, 3.2), sharex=True)
for a, g in zip(axs.ravel(), T.GROUP_ORDER):
    x = ana[ana.group == g].logA.to_numpy(); a.hist(x, bins=16, density=True, alpha=.7)
    t = np.linspace(x.min() - .3, x.max() + .3, 200); a.plot(t, stats.norm.pdf(t, x.mean(), x.std()), "C3", lw=1)
    a.set_title(f"{T.short_group(g)} (n={len(x)})", fontsize=7)
for a in axs[1]: a.set_xlabel(r"$\log A$")
fig.tight_layout(); fig.savefig(OUT / "logA_groups.pdf"); plt.close(fig)

# Fig 3: Delta with the three methods
fig, ax = plt.subplots(figsize=(3.4, 2.2))
f3, b = s3["nucleus_vs_cito"], s4["delta_posterior"]; h = s4["hierarchical"]
prof = {"MCF7R": (-0.8139, -0.8904, -0.7374), "T47DR": (-0.2468, -0.3425, -0.1510)}   # profile likelihood, notebook 03
for i, l in enumerate(lines):
    rows = [("profile likelihood", prof[l]), ("Bayesian", (b[l]["median"], b[l]["lo68"], b[l]["hi68"])),
            ("hierarchical", (h["b+u_" + l]["median"], h["b+u_" + l]["lo68"], h["b+u_" + l]["hi68"]))]
    for j, (nm, (m, lo, hi)) in enumerate(rows):
        y = i * 4 + j; ax.errorbar(m, y, xerr=[[m - lo], [hi - m]], fmt="o", color=f"C{j}", capsize=2, label=nm if i == 0 else None)
ax.set_yticks([1, 5], lines); ax.axvline(0, c="k", lw=.6); ax.set_xlabel(r"$\Delta=\mu_{nucleus}-\mu_{cytoplasm}$  (68%)"); ax.legend(fontsize=6, loc="center left")
fig.tight_layout(); fig.savefig(OUT / "delta_methods.pdf"); plt.close(fig)

# Fig 4: Bayes factor vs prior width and systematics
v = pd.read_csv(R / "fit_variants.csv"); v = v[np.isfinite(v.logA)]
fig, ax = plt.subplots(1, 2, figsize=(3.4, 2.2))
bf = pd.DataFrame(s4["log10_BF10"]); bf.index = bf.index.astype(float)
for l in lines: ax[0].plot(bf.index, bf[l], "o-", ms=3, label=l)
ax[0].axhline(0, c="k", lw=.6); ax[0].set_xscale("log"); ax[0].set_xlabel(r"prior width $\tau$"); ax[0].set_ylabel(r"$\log_{10}BF_{10}$"); ax[0].legend(fontsize=6)
def dlt(var, l):
    s = v[v.variant == var]; return (s[(s.cell_line if "cell_line" in s else s.group).astype(str).str.contains(l) & s.group.str.contains("nucleo")].logA.mean()
                                   - s[s.group.str.contains(l) & s.group.str.contains("cito")].logA.mean())
vars_ = [x for x in v.variant.unique() if x != "nominal"]
for i, l in enumerate(lines):
    ax[1].barh(np.arange(len(vars_)) + .4 * i, [dlt(x, l) - dlt("nominal", l) for x in vars_], .4, label=l)
ax[1].set_yticks(np.arange(len(vars_)) + .2, [x.replace("_", " ") for x in vars_], fontsize=5); ax[1].axvline(0, c="k", lw=.6); ax[1].set_xlabel(r"$\Delta_{var}-\Delta_{nom}$")
fig.tight_layout(); fig.savefig(OUT / "bf_systematics.pdf"); plt.close(fig)

# Fig 5: ML
fig, ax = plt.subplots(1, 2, figsize=(3.4, 2.2))
eff = pd.DataFrame(s5["anomaly"]["efficiency"]).T[["IsolationForest", "Robust Mahalanobis"]].astype(float)
eff.plot.bar(ax=ax[0], width=.8, legend=True); ax[0].set_ylabel("efficiency (5% FPR)"); ax[0].set_xticklabels([x.replace("_", "\n") for x in eff.index], rotation=0, fontsize=5)
ax[0].legend(fontsize=5)
cl = pd.DataFrame(s5["classification"]); piv = cl.pivot(index="line", columns="model", values="AUC")
piv.plot.bar(ax=ax[1], width=.8); ax[1].axhline(.5, c="k", lw=.6, ls=":"); ax[1].set_ylabel("AUC (CV)"); ax[1].set_xlabel(""); ax[1].legend(fontsize=5); ax[1].set_xticklabels(piv.index, rotation=0)
fig.tight_layout(); fig.savefig(OUT / "ml.pdf"); plt.close(fig)
print("figures written to", OUT)


# ---------------------------------------------------------------------------
# Extra figures: toy MC (pull, power), MCMC posterior checks, ROC/calibration.
# These re-run small versions of the notebook computations (a minute or two).
# ---------------------------------------------------------------------------
from scipy import optimize
import emcee
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.calibration import calibration_curve

X = {(l, r): ana[(ana.cell_line == l) & (ana.region == r)].logA.to_numpy() for l in lines for r in ("cito", "nucleo")}

def prof(delta, xc, xn):
    nc, nn = len(xc), len(xn)
    f = lambda mu: nc * np.log(np.mean((xc - mu) ** 2)) + nn * np.log(np.mean((xn - mu - delta) ** 2))
    lo, hi = min(xc.min(), (xn - delta).min()), max(xc.max(), (xn - delta).max())
    return optimize.minimize_scalar(f, bounds=(lo, hi), method="bounded", options={"xatol": 1e-9}).fun

def dfit(xc, xn):
    d = xn.mean() - xc.mean(); p0 = prof(d, xc, xn); f = lambda x: prof(x, xc, xn) - p0
    s = np.hypot(xc.std() / np.sqrt(len(xc)), xn.std() / np.sqrt(len(xn)))
    lo = optimize.brentq(lambda x: f(x) - 1, d - 8 * s, d); hi = optimize.brentq(lambda x: f(x) - 1, d, d + 8 * s)
    return d, (hi - lo) / 2, stats.chi2.sf(f(0.0), 1)

def toys(l, n=600, shift=None, seed=0):
    rg = np.random.default_rng(seed); xc0, xn0 = X[l, "cito"], X[l, "nucleo"]
    dt = xn0.mean() - xc0.mean() if shift is None else shift; out = []
    for _ in range(n):
        xc = rg.normal(xc0.mean(), xc0.std(), len(xc0)); xn = rg.normal(xc0.mean() + dt, xn0.std(), len(xn0)); out.append(dfit(xc, xn))
    return dt, np.array(out)

fig, ax = plt.subplots(1, 2, figsize=(3.4, 2.2)); shifts = np.linspace(0, .5, 9)
for l in lines:
    dt, a = toys(l); ax[0].hist((a[:, 0] - dt) / a[:, 1], bins=30, density=True, alpha=.5, label=l)
    ax[1].plot(shifts, [(toys(l, 250, s, 3)[1][:, 2] < .05).mean() for s in shifts], "o-", ms=3, label=l)
t = np.linspace(-4, 4, 200); ax[0].plot(t, stats.norm.pdf(t), "k", lw=.8); ax[0].set_xlabel(r"pull of $\Delta$"); ax[0].legend(fontsize=6)
ax[1].axhline(.8, c="gray", ls="--", lw=.6); ax[1].set_xlabel(r"true $|\Delta|$"); ax[1].set_ylabel("power (5%)")
fig.tight_layout(); fig.savefig(OUT / "toys.pdf"); plt.close(fig)

# MCMC (base model) -> trace, posterior predictive skewness
xs = [X[k] for k in [("MCF7R", "cito"), ("MCF7R", "nucleo"), ("T47DR", "cito"), ("T47DR", "nucleo")]]
def logp(th):
    mu, ls = th[:4], th[4:]; s = np.exp(ls)
    if np.any(ls < -6) or np.any(ls > 3): return -np.inf
    lp = -0.5 * np.sum(mu ** 2) + np.sum(-0.5 * s ** 2 + ls)
    return lp + sum(-0.5 * np.sum((x - m) ** 2) / si ** 2 - len(x) * np.log(si) for x, m, si in zip(xs, mu, s))
p0 = np.r_[[x.mean() for x in xs], np.log([x.std() for x in xs])]
rg = np.random.default_rng(1); sam = emcee.EnsembleSampler(40, 8, logp); sam.run_mcmc(p0 + 1e-2 * rg.standard_normal((40, 8)), 5000, progress=False)
ch = sam.get_chain(); flat = ch[1500:].reshape(-1, 8)
fig, ax = plt.subplots(1, 2, figsize=(3.4, 2.2))
ax[0].plot(ch[:, :5, 1], lw=.3); ax[0].set_xlabel("step"); ax[0].set_ylabel(r"$\mu_{MCF7R,nucleus}$")
names = ["MCF7R\ncyto", "MCF7R\nnucl", "T47DR\ncyto", "T47DR\nnucl"]; pv = []
for k, x in enumerate(xs):
    sk = np.array([stats.skew(rg.normal(flat[i, k], np.exp(flat[i, 4 + k]), len(x))) for i in rg.choice(len(flat), 1500)]); pv.append((sk > stats.skew(x)).mean())
ax[1].bar(range(4), pv, color=["C0", "C0", "C1", "C1"]); ax[1].axhline(.5, c="k", lw=.6); ax[1].axhline(.05, c="r", ls=":", lw=.6)
ax[1].set_xticks(range(4), names, fontsize=5); ax[1].set_ylabel("PPC p-value (skewness)")
fig.tight_layout(); fig.savefig(OUT / "mcmc_ppc.pdf"); plt.close(fig)

# ROC and calibration of the logistic classifier
MECH = ["logA", "max_indentation_um", "adhesion_nN", "adhesion_work_nN_um"]
data = fit[fit.fit_status == "ok"].dropna(subset=MECH); fig, ax = plt.subplots(1, 2, figsize=(3.4, 2.2))
for l in lines:
    d = data[data.cell_line == l]; y = (d.region == "nucleo").astype(int).to_numpy()
    p = cross_val_predict(make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000)), d[MECH], y, cv=StratifiedKFold(5, shuffle=True, random_state=0), method="predict_proba")[:, 1]
    fp, tp, _ = roc_curve(y, p); ax[0].plot(fp, tp, label=f"{l} ({roc_auc_score(y, p):.2f})")
    fr, mp = calibration_curve(y, p, n_bins=6, strategy="quantile"); ax[1].plot(mp, fr, "o-", ms=3, label=l)
for a in ax: a.plot([0, 1], [0, 1], "k:", lw=.6)
ax[0].set_xlabel("false positive rate"); ax[0].set_ylabel("efficiency"); ax[0].legend(fontsize=6)
ax[1].set_xlabel("predicted probability"); ax[1].set_ylabel("nucleus fraction"); ax[1].legend(fontsize=6)
fig.tight_layout(); fig.savefig(OUT / "roc_calibration.pdf"); plt.close(fig)
print("extra figures written")
