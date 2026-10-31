---
output:
  html_document: default
  pdf_document: default
---
# Supplementary Figures

## Fig. S — Pairwise distribution and correlation of the environmental covariates

![](figures/figures_manuscript/covariate_plot.pdf)

Pairwise matrix of the five environmental covariates used as predictors in the HMSC models (average NDVI, precipitation, maximum temperature, temperature range and wind speed), extracted at the 88 sampling points. Diagonal panels show the marginal distribution of each variable, lower panels show pairwise scatterplots, and upper panels report the corresponding Pearson correlation coefficients.

---

## Fig. S — AUC comparison between best- and second-best WAIC occurrence models

![](figures/figures_manuscript/auc_comparison_best_vs_second.pdf)

Per-genus cross-validated AUC for the two top occurrence models selected on WAIC (best in orange, second-best in blue). The dashed vertical line (AUC < 0.5) indicates predictive performance no better than chance.

---

## Fig. S — MCMC convergence diagnostics of the two retained HMSC models

![](figures/figures_manuscript/convergence_diagnostics.pdf)

Distribution of the Gelman–Rubin potential scale reduction factor (top row) and effective sample size (bottom row) of the Beta (species–environment) posterior, computed across 10 parallel chains for the retained abundance (Gaussian) and occurrence (probit) models. The dashed vertical line on the PSRF panels indicates the conventional 1.1 convergence threshold.

---

## Fig. S — DHARMa residual diagnostics of the retained abundance (Gaussian) model

![](figures/figures_manuscript/residuals_dharma_gauss_best_waic.svg)

Left: quantile–quantile plot of the scaled residuals against the uniform distribution, with Kolmogorov–Smirnov, dispersion and outlier tests. Right: scaled residuals against predicted values, with quantile regression lines (0.25, 0.5, 0.75) expected to align with the dashed reference lines under a correctly specified model.

---

## Fig. S — DHARMa residual diagnostics of the retained occurrence (probit) model

![](figures/figures_manuscript/residuals_dharma_bin_second_best_waic.svg)

Left: quantile–quantile plot of the scaled residuals against the uniform distribution, with Kolmogorov–Smirnov, dispersion and outlier tests. Right: scaled residuals against predicted values, with quantile regression lines (0.25, 0.5, 0.75) expected to align with the dashed reference lines under a correctly specified model.

---

## Fig. S — Spatial autocorrelation of the abundance-model residuals (DHARMa Moran's I)

![](figures/figures_manuscript/spatial_autocorr_dharma_gauss_second_best_waic.pdf)

Distance-based Moran's I test performed on site-aggregated DHARMa residuals of the retained Gaussian abundance model (observed I = −0.080, expected I = −0.011, sd = 0.081, p = 0.39). Blue and red dots indicate spatial coordinate with correlated residuals.

---

## Fig. S — Spatial autocorrelation of the occurrence-model residuals (DHARMa Moran's I)

![](figures/figures_manuscript/spatial_autocorr_dharma_bin_second_best_waic.pdf)

Distance-based Moran's I test performed on site-aggregated DHARMa residuals of the retained probit occurrence model (observed I = 0.064, expected I = −0.011, sd = 0.081, p = 0.35). Blue and red dots indicate spatial coordinate with correlated residuals.
