# Factorial Multivariate Bayesian Causal Forests — code and materials

Companion repository to the paper *Factorial Multivariate Bayesian Causal
Forests: heterogeneous main and interaction effects of multiple treatments on
correlated outcomes*. It contains everything needed to reproduce the results;
the editable manuscript is maintained as an Overleaf project (a compiled PDF is
included here for reference).

## Layout
- `mvbcfMT/` — the R package implementing the method (fit, per-component
  estimands, cross-fitted debiased average effects, overlap/sensitivity
  diagnostics, interpretability layer). C++ engine via Rcpp/RcppArmadillo.
- `code/` — reproducible analysis scripts: NHANES application
  (`nhanes_application.R`, `debias_nhanes.R`), the ACCORD template
  (`accord_application.R`), and `simulations/`.
- `results/` — average-effect tables (CSV), grf comparison, debiased estimates,
  figures.
- `paper/manuscript.pdf` — compiled reference copy (editable source on Overleaf).

## Reproducing
```r
R CMD INSTALL mvbcfMT
Rscript code/nhanes_application.R
Rscript code/debias_nhanes.R
```

## Built on
McJames, N., O'Shea, A., Goh, Y. C., & Parnell, A. (2025). Bayesian causal
forests for multivariate outcomes. JRSS Series A, 188(2), 428-450.
