# Mechanical stiffness of tumour cells from AFM force–distance curves

Analysis of AFM force–distance curves (JPK, MLCT-BIO cantilever) on two resistant breast-carcinoma cell lines
(MCF7R and T47DR), measured on the cytoplasm and on the nucleus at 48 h. University of Milano-Bicocca (Monza).

From the raw curves we extract the effective elastic parameter *A* of the Hertz-like model
$F(\delta)=A\,\delta^{3/2}$ and compare its distributions between regions and lines with frequentist, Bayesian and
machine-learning techniques, applying the methods of the *Statistics and Data Analysis* course
(same organisation as the `LHCb` project).

> *A* is an **effective** parameter in nN/µm^3/2, not the Young modulus: the tip radius and the Poisson ratio
> would be needed (absent from the JPK metadata). Comparisons are therefore relative.

Region labels follow the data folders: `cito` = cytoplasm, `nucleo` = nucleus.

## Repository layout

```
Code_General/
  01_JPK_reader_and_preliminary_QC.ipynb   read the .jpk-force files, export CSV, diagnostic instrumental QC
  02_preprocessing_contact_point.ipynb     baseline, contact point, Hertz-like fit, variants for the systematics
  03_frequentist_analysis.ipynb            MLE, profile likelihood, toy MC, interaction, Holm, systematics
  04_bayesian_analysis.ipynb               MCMC, Bayes factors, priors, PPC, WAIC, hierarchical, low statistics
  05_ml_curve_analysis.ipynb               anomaly detection (injected defects), classification, mixtures
  afm_tools.py                             shared functions (baseline, contact, fit, bootstrap, style)
  results/                                 light results, tracked: fit tables and summary JSON
  Preliminary/                             superseded exploratory prototypes (v0–v4), reference only
  Curve di forza belle/                    raw .jpk-force data (tracked, 78 MB)
  afm_output/                              cache regenerable from notebook 01 (NOT in the repo)
Docs/report/                               the report (main.tex, biblio.bib, figures, main.pdf); `make_figures.py` rebuilds the figures
Docs/                                      slides, QC guide, proof of concept on synthetic data
environment.yml                            conda environment of all the code
```

Only notebook 01 needs the raw data; 02 needs the cache `afm_output/` produced by 01.
Notebooks **03, 04 and 05 read only `results/`** (tracked) and run without the data.

## Data

The raw `.jpk-force` files (4 groups, 257 curves, 78 MB) are included in

```
Code_General/Curve di forza belle/<group>/*.jpk-force
```

with the groups `Petri20k_MCF7R_t48h_{cito,nucleo}_20260225` and `T47DR_t48h_{cito,nucleo}_20260318`. Running notebooks
01 → 02 in order regenerates the cache `afm_output/` (about 2.5 GB, not tracked) and `results/`. The folder name encodes cell line,
region, time and date.

## Environment

```
conda env create -f environment.yml
conda activate afm
python -m ipykernel install --user --name cern
```

The notebooks use the kernel `cern`. Approximate run times: 02 ≈ 2.5 min, 03 ≈ 2 min, 04 ≈ 3 min, 05 ≈ 1 min.

The complete report is `Docs/report/main.pdf` (build with `tectonic main.tex`).

## Pipeline

| Notebook | Content | Output |
|---|---|---|
| `01` | reads all `.jpk-force`, exports the segments to CSV, preliminary instrumental QC (robust z-scores). **No curve is discarded.** | `afm_output/` |
| `02` | robust linear baseline, contact point, fit $A\,\delta^{3/2}$ with free contact, quality flag, 9 preprocessing variants, adhesion from the retract | `results/fit_results_segment0.csv`, `results/fit_variants.csv` |
| `03` | unbinned MLE (`iminuit`), profile likelihood of Δ = μ_nucleus − μ_cytoplasm, comparison with Welch/Mann–Whitney/permutation/bootstrap, toy MC (bias, pull, coverage, size, power), interaction test, Holm, systematics (Barlow) | `results/03_summary.json` |
| `04` | posterior of Δ (`emcee`), Bayes factors versus prior width, prior sensitivity, posterior predictive check, Student-t (WAIC), hierarchical model, credible vs profile intervals at low statistics | `results/04_summary.json` |
| `05` | synthetic curves with 5 injected defects to measure the efficiency of Isolation Forest and robust Mahalanobis, cytoplasm/nucleus classification with CV and permutation, transfer between lines, Gaussian mixtures | `results/05_summary.json`, `results/05_real_anomaly_scores.csv` |

## Main results

Variable: log *A* (A in nN/µm^3/2); Δ = μ_nucleus − μ_cytoplasm, e^Δ = nucleus/cytoplasm ratio of the typical *A*.

| Line | n (cito/nucleo) | Δ (profile, 68%) | e^Δ | p (Wilks) | log₁₀ BF₁₀ (τ=0.5) | AUC cyto/nucleus |
|---|---|---|---|---|---|---|
| MCF7R | 55 / 60 | −0.81 (−0.89, −0.74) | 0.44 | < 10⁻¹⁵ | 13.9 | 0.95 |
| T47DR | 59 / 83 | −0.25 (−0.34, −0.15) | 0.78 | 0.011 | 0.6 | 0.63 |

- In both lines the nucleus is **softer** than the cytoplasm according to *A*; the effect is clear in MCF7R and
  weak in T47DR (the Bayes factor is ~4:1 despite p ≈ 0.01).
- The effect depends on the line (interaction test p ≈ 6·10⁻⁶; hierarchical τ_u ≈ 0.4).
- The preprocessing systematics (9 variants) stay ≤ 0.07 in Δ, comparable to the statistical error.
- The method is validated with toy MC: correct 68% coverage, size ≈ 5%, 80% power for |Δ| ≈ 0.22–0.25.
- The log *A* distribution of T47DR nucleus has a tail at high values (2-component mixture, p ≈ 0.01).

## Limits

- **Pseudo-replication:** several curves come from the same cell and the cell identifier is not available;
  the statistical errors are optimistic.
- **Line ≡ session:** MCF7R (25/02) and T47DR (18/03) were measured on different days: the comparison between
  lines cannot be separated from session effects. The nucleus–cytoplasm comparison within a line can.
- **Pre-selected sample:** the «good curves» all pass the QC; there is no sample of real defects.
  The anomaly-detection efficiency is estimated on synthetic curves and is optimistic.
- *A* is not the Young modulus and assumes the 3/2 exponent (sphere); the MLCT-BIO tip is pyramidal.

## Conventions

Notebooks and comments are in English. Every notebook states at the beginning what it does and what it does **not**
do, and ends with the output produced and the next notebook. `afm_output/` is a cache regenerable from notebook 01;
`results/` is the source of the results quoted here.
