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
  06_ml_shape_extensions.ipynb             curve-shape ML: PCA, clustering, conformal prediction, fixed-indentation check
  07_depth_resolved_stiffness.ipynb        stiffness vs indentation depth, fixed-force check, T47DR shape signal
  08_time_and_clustering_checks.ipynb      acquisition order, curve bursts (pseudo-replication), temporal drift
  09_calibration_and_repeatability.ipynb   calibration between sessions, repeatability, preconditioning
  10_unified_bayesian_model.ipynb          one hierarchical model: region, drift, position in the burst, burst effect
  afm_tools.py                             shared functions (baseline, contact, fit, bootstrap, style)
  results/                                 light results, tracked: fit tables and summary JSON
  Preliminary/                             superseded exploratory prototypes (v0–v4), reference only
  Curve di forza belle/                    raw .jpk-force data (NOT in the repo, see below)
  afm_output/                              cache regenerable from notebook 01 (NOT in the repo)
Docs/report/                               the report (main.tex, biblio.bib, figures, main.pdf); `make_figures.py` rebuilds the figures
environment.yml                            conda environment of all the code
```

Only notebook 01 needs the raw data; 02 needs the cache `afm_output/` produced by 01.
Notebooks **03–10 read only `results/`** (tracked) and run without the data.

## Data

The raw `.jpk-force` files (4 groups, 257 curves) belong to the laboratory that measured them and are **not** distributed with
the repository. To regenerate everything they must be placed in

```
Code_General/Curve di forza belle/<group>/*.jpk-force
```

with the groups `Petri20k_MCF7R_t48h_{cito,nucleo}_20260225` and `T47DR_t48h_{cito,nucleo}_20260318`, then notebooks 01 → 02
are run in order. The folder name encodes cell line, region, time and date. Notebooks 03–05 need only `results/`, which is included.

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
| `06` | curve-shape extensions: PCA of the force–indentation curve, clustering, classification from shape vs scalar features, force at fixed indentation as model-independent check, conformal prediction, learning curves, shape-based anomaly detector | `results/06_summary.json`, `results/curves_resampled.csv`, `results/06_real_shape_anomaly.csv` |
| `07` | depth-resolved nucleus–cytoplasm difference (force at fixed indentation, indentation at fixed force, effective stiffness vs depth) and test of the T47DR shape signal against contact-point shifts | `results/07_summary.json`, `results/force_indentation_profiles.csv` |
| `08` | acquisition timestamps: order of measurement, bursts of consecutive curves as a proxy of the cell (ICC, design effect), burst-level tests, temporal drift and its effect on the region difference | `results/08_summary.json`, `results/acquisition_times.csv` |
| `09` | calibration metadata per session (the stored spring constant of the same cantilever differs by 40%), variance components at the same spot, preconditioning (first curve at a spot is stiffer), bias of the region difference from burst length | `results/09_summary.json`, `results/session_calibration.csv` |
| `10` | one Bayesian hierarchical model with region, time drift, position in the burst and burst random effect, for log A, force at 300 nm and free exponent, with prior sensitivity of the drift | `results/10_summary.json` |

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
- Model-independent checks (force at fixed indentation, notebooks 06–07) give the same sign in both lines, but see the caveats below: for T47DR they are subject to the same temporal confounding as *A*. The free contact-law exponent is ≈ 1.5, consistent with the model used. At the level of the raw comparison the nucleus is softer at every indentation depth and needs a larger indentation to reach the same force (notebook 07). In T47DR the nuclear curves also have a higher free exponent (1.50 vs 1.43, p ≈ 5·10⁻⁸): the region differs in the form of the response, not only in stiffness.

## Caveats from the acquisition times and calibration (notebooks 08–10)

- In both sessions the cytoplasm was measured first and the nucleus afterwards, in non-overlapping time blocks: region is confounded with time.
- The stored spring constant of the same cantilever differs by 40% between the two sessions (0.119 vs 0.080 N/m), so differences **between lines** cannot be interpreted, not even in sign. Within a session the calibration cancels in the nucleus–cytoplasm difference.
- The first curve at a spot is stiffer than the following ones (second − first = −0.25 ± 0.09 in log *A*), and the nucleus bursts are longer than the cytoplasm ones.
- MCF7R: the nucleus is softer than the cytoplasm in every model (β ≈ −0.8 to −0.9 in log *A*, 95% interval [−1.2, −0.55]), robust to drift, bursts and preconditioning.
- T47DR: the stiffness difference is **not established**: β = −0.22 [−0.45, +0.02] with bursts only, +0.03 [−0.35, +0.41] with the drift term, and ≈ 0 using only the first curve of each burst. The higher free exponent of the T47DR nucleus (+0.06) is robust to all these corrections.

## Limits

- **Pseudo-replication:** several curves come from the same cell and the cell identifier is not available;
  the statistical errors are optimistic.
- **Line ≡ session:** MCF7R (25/02) and T47DR (18/03) were measured on different days: the comparison between
  lines cannot be separated from session effects. The nucleus–cytoplasm comparison within a line can.
- **Pre-selected sample:** the «good curves» all pass the QC; there is no sample of real defects.
  The anomaly-detection efficiency is estimated on synthetic curves and is optimistic.
- *A* is not the Young modulus and assumes the 3/2 exponent (sphere); the free exponent estimated from the data is ≈ 1.5, but the tip geometry is not in the metadata.

## Conventions

Notebooks and comments are in English. Every notebook states at the beginning what it does and what it does **not**
do, and ends with the output produced and the next notebook. `afm_output/` is a cache regenerable from notebook 01;
`results/` is the source of the results quoted here.

## License

The code and the report are released under the MIT License (see `LICENSE`). The raw JPK force curves are
**not** covered by this license: they belong to the laboratory that measured them and are not distributed here.
