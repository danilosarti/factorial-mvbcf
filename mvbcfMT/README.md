# mvbcfMT

<!-- badges: start -->
<!-- badges: end -->

**Multi-Treatment Multivariate Bayesian Causal Forests.**

`mvbcfMT` estimates the joint causal effect of two or more crossed binary
treatments on several correlated continuous outcomes at once. It extends the
single-treatment multivariate Bayesian causal forest to the full factorial: the
response surface is decomposed over the treatment lattice into a prognostic
sum-of-trees forest plus one forest per main or interaction effect, all sharing
an inverse-Wishart residual covariance across outcomes. A single unified `Rcpp`
engine samples every component, so the two-treatment model, the single-treatment
multivariate causal forest, and a general order-*r* truncation for more
treatments are one and the same sampler.

The package provides:

- per-treatment propensity adjustment in the prognostic forest;
- per-component main, interaction and joint effect estimands with credible
  intervals (`component_effect()`, `contrast_effect()`, `regime_outcome()`,
  `predict_effects()`);
- a cross-fitted debiased (one-step / AIPW) average-effect estimator with
  efficient standard errors and doubly-robust Wald intervals (`debiased_ate()`);
- positivity / factorial-overlap and unmeasured-confounding diagnostics
  (`overlap_diagnostic()`, `sensitivity_effect()`);
- value-suppressing uncertainty maps of the interaction surface
  (`plot_effect_vsup_xy()`, `plot_effect_forest()`, `vivi_matrix()`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("daniloasarti/mvbcfMT")  # once the repository is public
```

The package compiles C++ via `Rcpp`, `RcppArmadillo` and `RcppDist`; a working
C++17 toolchain is required.

## Quick start

```r
library(mvbcfMT)

## simulate two crossed treatments, two correlated outcomes
sim <- simulate_multi(n = 300, q = 2, Ti = 2, seed = 1)

## fit the factorial multivariate BCF
fit <- fit_mvbcf_multi(
  sim$data,
  responses  = sim$responses,
  treatments = sim$treatments,
  covariates = sim$covariates,
  n_iter = 500, n_burn = 250
)

## average main / interaction / joint effects, calibrated
predict_effects(fit, test = FALSE)

## positivity and overlap before trusting any interaction
overlap_diagnostic(fit)

## semiparametric-efficient, doubly-robust average effects
debiased_ate(fit)
```

## Discipline

The population **average** effect (ATE) is the calibrated primary claim;
per-unit (CATE) surfaces are exploratory and shown with value-suppressing
uncertainty maps. An interaction is reported only where its factorial cell has
support — the overlap diagnostic flags thin cells as *unsupported*, not *null*.

## Citation

See `citation("mvbcfMT")`.

## License

MIT © Danilo A. Sarti.
