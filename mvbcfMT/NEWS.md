# mvbcfMT 0.1.0

* First release.
* Factorial multivariate Bayesian causal forest for two or more crossed binary
  treatments on several correlated continuous outcomes, sampled by a single
  unified `Rcpp` engine (`fit_mvbcf_multi()`, `fit_mvbcf2()`).
* Per-component estimands with credible intervals: `component_effect()`,
  `contrast_effect()`, `regime_outcome()`, `predict_effects()`.
* Cross-fitted debiased (one-step / AIPW) average-effect estimator with
  efficient standard errors and doubly-robust Wald intervals: `debiased_ate()`.
* Positivity / factorial-overlap and unmeasured-confounding diagnostics:
  `overlap_diagnostic()`, `sensitivity_effect()`, `calibrate_intervals()`.
* Value-suppressing uncertainty interpretability layer: `plot_effect_vsup_xy()`,
  `plot_effect_forest()`, `vivi_matrix()`, `plot_vivi_vsup()`.
* Data generators for validation: `simulate_multi()`, `simulate_met2()`.
