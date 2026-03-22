# Nonparametric Spatio-Temporal Hawkes Modeling of Lightning Events

This repository contains R code for detecting lightning-event clusters, fitting a nonparametric spatio-temporal Hawkes model, and evaluating model fit with super-thinning residuals and related summaries.

The codebase is organized around a simple workflow:

1. detect lightning-event windows from raw lightning data,
2. compare or visualize alternative event definitions,
3. fit a nonparametric spatio-temporal Hawkes model to each event window,
4. compute residual diagnostics and summary plots.

---

## Repository structure

```text
st-hawkes-lightning/
├─ R/
│  ├─ preprocessing/
│  │  └─ detect_clusters_by_threshold.R
│  ├─ EDA/
│  │  ├─ compare_thresholds.R
│  │  ├─ run_stdbscan.R
│  │  └─ summarize_clusters.R
│  ├─ methods/
│  │  ├─ fit_princurve.R
│  │  ├─ run_hawkes_by_year.R
│  │  ├─ run_hawkes_year_threshold_grid.R
│  │  ├─ super_thinning.R
│  │  └─ summaries_and_plots.R
├─ plots/
└─ README.md
