#' Posterior summaries (mean + CI) of one component on train or test design.
#' @keywords internal
.component_post <- function(fit, k, test = TRUE, probs = c(0.025, 0.975)) {
  cube <- if (test) fit$model$predictions_test[[k]] else fit$model$predictions[[k]]
  list(mean = apply(cube, c(1, 2), mean),
       lo   = apply(cube, c(1, 2), stats::quantile, probs[1]),
       hi   = apply(cube, c(1, 2), stats::quantile, probs[2]),
       draws = cube)
}

#' Extract a component's posterior effect surface
#'
#' @param fit A \code{\link{fit_mvbcf_multi}} object.
#' @param component Component label (e.g. \code{"tau_Z1"}, \code{"tau_Z1:Z2"},
#'   \code{"mu"}) or index.
#' @param test Use the test/grid design (default) or the training rows.
#' @return A list with \code{mean}, \code{lo}, \code{hi} matrices (rows x q).
#' @examples
#' \donttest{
#' sim <- simulate_multi(n = 150, q = 2, Ti = 2, seed = 1)
#' fit <- fit_mvbcf_multi(sim$data, responses = sim$responses,
#'                        treatments = sim$treatments, covariates = sim$covariates,
#'                        n_iter = 100, n_burn = 50, n_tree = 30, n_tree_tau = 20)
#' eff <- component_effect(fit, "tau_Z1", test = FALSE)
#' colMeans(eff$mean)
#' }
#' @export
component_effect <- function(fit, component, test = TRUE) {
  k <- if (is.numeric(component)) component else match(component, fit$components$labels)
  if (is.na(k)) stop("unknown component; available: ",
                     paste(fit$components$labels, collapse = ", "))
  .component_post(fit, k, test = test)
}

#' Posterior draws of an arbitrary treatment-configuration contrast
#'
#' Computes, per row and per posterior draw, the potential-outcome contrast
#' between two treatment configurations \code{a} and \code{b} (each a 0/1 vector
#' over the treatments). The engine stores every \eqn{\tau_S} surface, so any
#' contrast is a signed sum of components:
#' \deqn{Y(a)-Y(b) = \sum_{\emptyset\ne S} \Big(\prod_{t\in S}a_t - \prod_{t\in S}b_t\Big)\tau_S.}
#' Main effect of treatment t = contrast(e_t, 0); interaction of (s,t) =
#' contrast with the 2x2 difference-in-differences weights; joint effect =
#' contrast(1, 0).
#'
#' @param fit A fit object.
#' @param a,b Named or positional 0/1 vectors over \code{fit$treatments}.
#' @param test Use the test/grid design (default) or training rows.
#' @return A list: \code{mean}, \code{lo}, \code{hi} (rows x q) and \code{draws}.
#' @export
contrast_effect <- function(fit, a, b = rep(0, length(fit$treatments)), test = TRUE) {
  Ti <- length(fit$treatments); a <- as.numeric(a); b <- as.numeric(b)
  subs <- fit$components$subsets
  first <- if (test) fit$model$predictions_test[[1]] else fit$model$predictions[[1]]
  acc <- array(0, dim = dim(first))
  for (k in seq_along(subs)) {
    S <- subs[[k]]; if (length(S) == 0) next
    w <- prod(a[S]) - prod(b[S]); if (w == 0) next
    cube <- if (test) fit$model$predictions_test[[k]] else fit$model$predictions[[k]]
    acc <- acc + w * cube
  }
  list(mean = apply(acc, c(1, 2), mean),
       lo   = apply(acc, c(1, 2), stats::quantile, 0.025),
       hi   = apply(acc, c(1, 2), stats::quantile, 0.975),
       draws = acc)
}

#' Predicted outcome level under a given treatment configuration
#'
#' \eqn{E[Y(z)] = \mu + \sum_{S \subseteq \mathrm{on}(z)} \tau_S}.
#' @param fit A fit object.
#' @param z 0/1 vector over the treatments.
#' @param test test/grid design (default) or training rows.
#' @return list of \code{mean}, \code{lo}, \code{hi}, \code{draws}.
#' @export
regime_outcome <- function(fit, z, test = TRUE) {
  z <- as.numeric(z); subs <- fit$components$subsets
  first <- if (test) fit$model$predictions_test[[1]] else fit$model$predictions[[1]]
  acc <- array(0, dim = dim(first))
  for (k in seq_along(subs)) {
    S <- subs[[k]]
    w <- if (length(S) == 0) 1 else prod(z[S])   # mu always on; tau_S on iff all t in S active
    if (w == 0) next
    cube <- if (test) fit$model$predictions_test[[k]] else fit$model$predictions[[k]]
    acc <- acc + w * cube
  }
  list(mean = apply(acc, c(1, 2), mean),
       lo   = apply(acc, c(1, 2), stats::quantile, 0.025),
       hi   = apply(acc, c(1, 2), stats::quantile, 0.975),
       draws = acc)
}

#' Tidy table of the main causal estimands for every unit / grid cell
#'
#' For two treatments returns, per response: the two main effects
#' (\eqn{\tau_1,\tau_2}), the interaction (\eqn{\tau_{12}}), the joint effect
#' (both on vs both off) and the outcome levels under all four regimes, each with
#' 95\% credible intervals. For more treatments returns \eqn{\mu}, each main
#' effect and (if fitted) each interaction component.
#'
#' @param fit A fit object.
#' @param test test/grid (default) or training rows.
#' @return A data frame, one row per unit/cell.
#' @examples
#' \donttest{
#' sim <- simulate_multi(n = 150, q = 2, Ti = 2, seed = 1)
#' fit <- fit_mvbcf_multi(sim$data, responses = sim$responses,
#'                        treatments = sim$treatments, covariates = sim$covariates,
#'                        n_iter = 100, n_burn = 50, n_tree = 30, n_tree_tau = 20)
#' head(predict_effects(fit, test = FALSE))
#' }
#' @export
predict_effects <- function(fit, test = TRUE) {
  resp <- fit$responses; q <- length(resp); tr <- fit$treatments
  base <- if (!is.null(fit$grid) && test)
    fit$grid[, c(fit$gen, fit$env, "observed")] else
    data.frame(row_id = seq_len(if (test) dim(fit$model$predictions_test[[1]])[1]
                                else nrow(fit$data)))
  out <- data.frame(base, stringsAsFactors = FALSE)
  add <- function(pref, e) for (k in seq_len(q)) {
    out[[paste0(pref, "_", resp[k])]]      <<- e$mean[, k]
    out[[paste0(pref, "_", resp[k], "_lo")]] <<- e$lo[, k]
    out[[paste0(pref, "_", resp[k], "_hi")]] <<- e$hi[, k]
  }
  # every fitted component surface
  for (lab in fit$components$labels) add(lab, component_effect(fit, lab, test = test))
  # joint (all-on vs all-off) for convenience
  add("joint_all", contrast_effect(fit, rep(1, length(tr)), rep(0, length(tr)), test = test))
  out
}
