#' mvbcfMT: Multi-Treatment Multivariate Bayesian Causal Forests
#'
#' Extends the single-treatment multivariate BCF (mvbcfMET) to estimate the joint
#' causal effect of \strong{two or more} crossed binary treatments on several
#' outcomes at once. The response surface is decomposed over the treatment
#' lattice into a prognostic forest \eqn{\mu} plus one forest per (main or
#' interaction) effect \eqn{\tau_S}, all sharing a residual covariance \eqn{\Sigma}.
#' A single unified C++ engine (\code{fast_bart_multi}) samples every component.
#'
#' @useDynLib mvbcfMT, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"
