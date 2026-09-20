## Submission summary

This is the first submission of mvbcfMT (0.1.0).

The package implements the factorial multivariate Bayesian causal forest that
accompanies the methods paper "Factorial Multivariate Bayesian Causal Forests"
(Sarti, 2026). It estimates heterogeneous main and interaction effects of several
crossed binary treatments on several correlated continuous outcomes, with a
single unified compiled sampler.

## Test environments

<!-- Fill in after running the checks locally: -->
* local: <your OS>, R <version>
* win-builder: devel and release (devtools::check_win_devel(), check_win_release())
* R-hub: windows, macos, ubuntu (rhub::rhub_check())
* macOS builder: macbuilder (devtools::check_mac_release())

## R CMD check results

<!-- Paste the summary of `R CMD check --as-cran` here. Target: -->
0 errors | 0 warnings | 0 notes

Expected/acceptable notes on a first submission:

* "New submission" (this is a new package).
* Possible timing note on the compiled examples; the average-effect and
  debiasing examples use tiny MCMC settings and are additionally wrapped in
  \donttest{}.

## Compiled code

* C++ via Rcpp, RcppArmadillo and RcppDist; standard set in src/Makevars
  (CXX_STD = CXX17).
* All console output from the sampler uses Rcpp::Rcout (no printf / std::cout).
* Random numbers are drawn through R's RNG (R::runif), so results respect
  set.seed() and no separate RNG stream is used. No OpenMP parallel regions call
  the R API.

## Downstream dependencies

None (first submission).
