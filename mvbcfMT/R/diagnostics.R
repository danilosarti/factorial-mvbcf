# Diagnostics and improvements motivated by the robustness study:
#  - overlap_diagnostic():   positivity / factorial-cell support (interaction id.)
#  - sensitivity_effect():   how strong an UNMEASURED confounder must be to nullify
#  - calibrate_intervals():  heuristic CATE-interval widening (undercoverage fix)

#' Overlap / positivity diagnostic
#'
#' The robustness study showed main effects survive strong \emph{measured}
#' confounding but interactions need every factorial cell populated. This reports,
#' per treatment, the estimated-propensity range and the share of near-violations
#' (p < 0.05 or > 0.95); and, per fitted interaction subset, how many units sit in
#' its "all-on" cell (\eqn{\prod_{t\in S}Z_t=1}) --- interactions whose on-cell is
#' thin are flagged as prior-driven, not data-driven.
#'
#' @param fit A \code{\link{fit_mvbcf_multi}} object.
#' @param thin Minimum on-cell count below which an interaction is flagged.
#' @return A list with \code{treatments} and \code{interactions} data frames.
#' @examples
#' \donttest{
#' sim <- simulate_multi(n = 200, q = 2, Ti = 2, seed = 1)
#' fit <- fit_mvbcf_multi(sim$data, responses = sim$responses,
#'                        treatments = sim$treatments, covariates = sim$covariates,
#'                        n_iter = 100, n_burn = 50, n_tree = 30, n_tree_tau = 20)
#' overlap_diagnostic(fit)
#' }
#' @export
overlap_diagnostic <- function(fit, thin = 30) {
  ps <- fit$propensity; tr <- fit$treatments
  Zmat <- as.matrix(fit$data[, tr, drop = FALSE])
  tdf <- data.frame(
    treatment = tr,
    p_min = round(apply(ps, 2, min), 3), p_max = round(apply(ps, 2, max), 3),
    frac_near_violation = round(colMeans(ps < 0.05 | ps > 0.95), 3),
    n_treated = colSums(Zmat == 1))
  subs <- fit$components$subsets; labs <- fit$components$labels
  irows <- lapply(seq_along(subs), function(k) {
    S <- subs[[k]]; if (length(S) < 2) return(NULL)
    oncell <- sum(apply(Zmat[, S, drop = FALSE], 1, prod) == 1)
    data.frame(interaction = labs[k], on_cell_n = oncell,
               supported = oncell >= thin)
  })
  idf <- do.call(rbind, irows)
  structure(list(treatments = tdf, interactions = idf), class = "mvbcf_overlap")
}

#' @export
print.mvbcf_overlap <- function(x, ...) {
  cat("Overlap / positivity diagnostic\n\n Treatments:\n"); print(x$treatments, row.names = FALSE)
  if (!is.null(x$interactions)) { cat("\n Interaction support (all-on cell):\n")
    print(x$interactions, row.names = FALSE)
    if (any(!x$interactions$supported))
      cat("  ! flagged interactions are prior-driven (thin cell) -- report as unsupported, not null.\n") }
  invisible(x)
}

#' Sensitivity of an effect to an unmeasured confounder
#'
#' No causal method survives unmeasured confounding (the study's \code{unmeasured}
#' scenario collapsed to ~0 coverage, as it must). This quantifies fragility: the
#' spurious effect an omitted confounder would have to induce to explain the point
#' estimate away, and (if the interval excludes zero) to reach the null boundary,
#' expressed in outcome units and in outcome-SD units.
#'
#' @param estimate,lo,hi Effect estimate and 95\% interval (same units).
#' @param outcome_sd SD of the outcome (to express bias in SD units).
#' @return A one-row data frame; \code{robustness_sd} is the omitted-confounder
#'   induced shift (in outcome SDs) needed to reach the nearest null boundary.
#' @export
sensitivity_effect <- function(estimate, lo, hi, outcome_sd = 1) {
  crosses <- lo <= 0 & hi >= 0
  bias_to_point <- abs(estimate)
  bias_to_null  <- if (crosses) 0 else min(abs(lo), abs(hi))
  data.frame(estimate = estimate, ci_lo = lo, ci_hi = hi,
             excludes_zero = !crosses,
             bias_to_nullify_point = round(bias_to_point, 4),
             bias_to_reach_CI_null = round(bias_to_null, 4),
             robustness_sd = round(bias_to_null / outcome_sd, 3))
}

#' Heuristic calibration of CATE credible intervals
#'
#' The study found per-unit (CATE) intervals under-cover (~0.70-0.90 vs 0.95),
#' while ATE intervals are calibrated. This widens the per-unit intervals of a
#' component by a global factor \eqn{\gamma} chosen so the model's \emph{outcome}
#' predictive intervals cover the observed responses at the nominal level --- a
#' data-driven (truth-free) inflation. It is a pragmatic partial fix, not a
#' guarantee; more trees/iterations and a heavier-tailed leaf prior address the
#' root cause (see the improvements note).
#'
#' @param fit A fit object.
#' @param component Component label/index.
#' @param response Which response.
#' @param level Nominal level.
#' @return A list: \code{gamma} and calibrated \code{lo}/\code{hi} (rows x 1).
#' @export
calibrate_intervals <- function(fit, component, response = fit$responses[1], level = 0.95) {
  ri <- match(response, fit$responses)
  # predictive spread of the fitted mean for the outcome (mu + active taus on train)
  z <- stats::qnorm(1 - (1 - level) / 2)
  # residual sd from stored Sigma draws
  sig <- sqrt(mean(sapply(seq_len(dim(fit$model$sigmas)[3]),
                          function(s) fit$model$sigmas[ri, ri, s])))
  # fitted mean per unit (sum of active components) posterior mean + its sd
  n <- nrow(fit$data); fitmean <- rep(0, n); fitvar <- rep(0, n)
  Zmat <- as.matrix(fit$data[, fit$treatments, drop = FALSE])
  for (k in seq_along(fit$components$subsets)) {
    S <- fit$components$subsets[[k]]
    d <- if (length(S) == 0) rep(1, n) else apply(Zmat[, S, drop = FALSE], 1, prod)
    cube <- fit$model$predictions[[k]][, ri, ]
    fitmean <- fitmean + d * apply(cube, 1, mean)
    fitvar  <- fitvar  + (d^2) * apply(cube, 1, stats::var)
  }
  y <- fit$data[[response]]
  # find gamma so predictive coverage of y hits `level`
  cover <- function(g) mean(abs(y - fitmean) <= z * sqrt(g^2 * fitvar + sig^2))
  gseq <- seq(1, 6, by = 0.1); g <- gseq[which.min(abs(sapply(gseq, cover) - level))]
  k <- if (is.numeric(component)) component else match(component, fit$components$labels)
  cube <- fit$model$predictions[[k]][, ri, ]
  m <- apply(cube, 1, mean); s <- apply(cube, 1, stats::sd)
  list(gamma = g, lo = m - z * g * s, hi = m + z * g * s,
       note = "CATE intervals widened by gamma (from outcome predictive calibration)")
}
