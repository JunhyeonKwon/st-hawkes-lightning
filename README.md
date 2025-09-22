## Nonparametric Spatio-Temporal Hawkes (Lightning)

This repo contains code to fit a **nonparametric spatio-temporal Hawkes process** (ETAS-type) to lightning events, compute **super-thinning residuals**, and evaluate **goodness-of-fit** with the spatio-temporal paired correlation function (PCF).

The code is designed to be **reproducible and modular**:

* `fit_nonpar_sthawkes.R` fits the model for each dataset window via an EM-style routine with step-function kernels.
* `super_thinning.R` computes super-thinning residuals and the residual PCF.
* `summaries_and_plots.R` makes publication-ready figures and quick diagnostics.
