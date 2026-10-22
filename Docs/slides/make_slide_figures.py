"""Figures specific to the slides (timeline, preconditioning, models). Run: python make_slide_figures.py"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
CODE = Path(__file__).resolve().parents[2] / "Code_General"; sys.path.insert(0, str(CODE)); import afm_tools as T
OUT = Path(__file__).resolve().parent / "figures"; T.set_plot_style(); plt.rcParams.update({"font.size": 10, "figure.dpi": 150})
R = CODE / "results"; lines = ["MCF7R", "T47DR"]
fit = T.load_fit_table(); fit = fit[fit.fit_status == "ok"]
times = pd.read_csv(R / "acquisition_times.csv", parse_dates=["t"]); df = fit.merge(times, on=["file_name", "group"]).sort_values(["cell_line", "region", "t"]).reset_index(drop=True)
df["t_min"] = df.groupby("cell_line").t.transform(lambda s: (s - s.min()).dt.total_seconds() / 60)
df["gap"] = df.groupby(["cell_line", "region"]).t.diff().dt.total_seconds(); df["burst"] = df.groupby(["cell_line", "region"], group_keys=False).apply(lambda d: ((d.gap.isna()) | (d.gap > 10)).cumsum())
df["cl"] = df.cell_line + df.region + df.burst.astype(str); df["pos"] = df.groupby("cl").cumcount() + 1

# timeline
fig, ax = plt.subplots(1, 2, figsize=(7.5, 3.6), sharey=True)
for a, l in zip(ax, lines):
    d = df[df.cell_line == l]
    for r_, c, lab in [("cito", "C0", "cytoplasm"), ("nucleo", "C1", "nucleus")]:
        s = d[d.region == r_]; a.scatter(s.t_min, s.logA, s=10, c=c, alpha=.8, label=lab); b = np.polyfit(s.t_min, s.logA, 1); tt = np.array([s.t_min.min(), s.t_min.max()]); a.plot(tt, np.polyval(b, tt), c=c)
    a.set_title(l); a.set_xlabel("minutes from the start of the session")
ax[0].set_ylabel(r"$\log A$"); ax[0].legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "timeline.pdf"); plt.close(fig)

# preconditioning
big = df[df.groupby("cl").logA.transform("size") >= 3].copy(); big["dev"] = big.logA - big.groupby("cl").logA.transform("mean"); rg = np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(4.2, 3.2)); ps = [1, 2, 3, 4]; m, lo, hi = [], [], []
cl = big.cl.unique()
for p in ps:
    m.append(big[big.pos == p].dev.mean()); bs = [big[(big.pos == p) & big.cl.isin(rg.choice(cl, len(cl)))].dev.mean() for _ in range(500)]; lo.append(np.nanquantile(bs, .025)); hi.append(np.nanquantile(bs, .975))
m, lo, hi = map(np.array, (m, lo, hi)); ax.errorbar(ps, m, [m - lo, hi - m], marker="o", capsize=3); ax.axhline(0, c="k", lw=.7)
ax.set_xlabel("position in the burst"); ax.set_ylabel(r"$\log A$ − burst mean"); ax.set_xticks(ps); fig.tight_layout(); fig.savefig(OUT / "preconditioning.pdf"); plt.close(fig)

# beta from the unified models
s10 = json.load(open(R / "10_summary.json")); tab = pd.DataFrame(s10["models"]); tab = tab[tab.observable == "logA"]
fig, ax = plt.subplots(figsize=(5.2, 3.8))
for j, l in enumerate(lines):
    t = tab[tab.line == l]
    for i, (_, r_) in enumerate(t.iterrows()):
        y = j * 5 + i; ax.errorbar(r_.beta_median, y, xerr=[[r_.beta_median - r_.lo95], [r_.hi95 - r_.beta_median]], fmt="o", color=f"C{i}", capsize=2, label=r_.model if j == 0 else None)
ax.axvline(0, c="k", lw=.7); ax.set_yticks([1.5, 6.5], lines); ax.set_xlabel(r"$\beta$ (nucleus − cytoplasm in $\log A$, 95%)"); ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2); fig.tight_layout(); fig.savefig(OUT / "beta_models.pdf"); plt.close(fig)

# drift-prior sensitivity
sens = pd.DataFrame(s10["drift_prior_sensitivity"]); fig, ax = plt.subplots(figsize=(4.4, 3.2))
for j, l in enumerate(lines):
    t = sens[sens.line == l]; ax.errorbar(t.s_gamma_sd_per_10min * (1 + .05 * j), t.beta_median, [t.beta_median - t.lo95, t.hi95 - t.beta_median], marker="o", capsize=2, label=l)
ax.axhline(0, c="k", lw=.7); ax.set_xscale("log"); ax.set_xlabel(r"prior width on the drift (sd per 10 min)"); ax.set_ylabel(r"$\beta$ for $\log A$ (95%)"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "drift_sensitivity.pdf"); plt.close(fig)
print("slide figures written")
