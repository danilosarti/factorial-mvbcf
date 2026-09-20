#' Cross-fitted debiased (one-step / AIPW) average effects for factorial MVBCF
#'
#' Implements the debiased average-effect estimator of the accompanying paper:
#' it corrects the MVBCF posterior-mean regime surfaces with a cross-fitted
#' augmented-IPW term over the 2^T factorial cells, yielding \eqn{\sqrt n}-efficient,
#' doubly-robust average main/interaction/joint effects with efficient standard
#' errors and Wald intervals -- discharging the no-bias condition of the semiparametric
#' Bernstein--von Mises result rather than assuming it.
#'
#' Nuisances: the outcome surfaces \eqn{\hat g_{S'}(x)=E[Y(1_{S'})\mid x]} come from
#' the fitted model via \code{regime_outcome()}; the cell propensities
#' \eqn{\hat e_{S'}(x)} are cross-fitted (ranger > nnet::multinom > per-cell glm).
#' Validity needs the product rate \eqn{\|\hat g-g\|\,\|\hat e-e\|=o_P(n^{-1/2})}
#' (cross-fitting removes the Donsker condition; see Chernozhukov et al. 2018,
#' Kennedy 2023).
#'
#' @param fit A \code{fit_mvbcf_multi} object with exactly two treatments (the 2x2).
#'   For general T the same code applies with the full Mobius combination of cells.
#' @param covariates Covariate columns for the propensity model; if \code{NULL},
#'   inferred from the fit's data (numeric columns, minus responses/treatments/ids).
#' @param folds Cross-fitting folds (default 5).
#' @param propensity One of "auto","ranger","multinom","glm".
#' @param trim Propensity trimming bound (default 0.01) for positivity.
#' @param level Credible/confidence level (default 0.95).
#' @param seed RNG seed for the fold split and forest propensity.
#' @return A data frame: for each component (tau_Z1, tau_Z2, tau_Z1:Z2, joint) and
#'   outcome, the raw posterior-plug-in ATE, the debiased estimate, its efficient SE,
#'   the Wald interval and an excludes-zero flag.
#' @examples
#' \donttest{
#' sim <- simulate_multi(n = 300, q = 2, Ti = 2, seed = 1)
#' fit <- fit_mvbcf_multi(sim$data, responses = sim$responses,
#'                        treatments = sim$treatments, covariates = sim$covariates,
#'                        n_iter = 100, n_burn = 50, n_tree = 30, n_tree_tau = 20)
#' debiased_ate(fit, propensity = "glm", folds = 3)
#' }
#' @export
debiased_ate <- function(fit, covariates = NULL, folds = 5,
                         propensity = c("auto","ranger","multinom","glm"),
                         trim = 0.01, level = 0.95, seed = 1) {
  propensity <- match.arg(propensity)
  tr <- fit$treatments; resp <- fit$responses; q <- length(resp)
  stopifnot(length(tr) == 2L)                       # 2x2 implementation
  d <- fit$data; n <- nrow(d)
  Z <- as.matrix(d[, tr]); storage.mode(Z) <- "integer"
  if (is.null(covariates)) {
    drop <- c(resp, tr, grep("^(pro_|true_)", names(d), value = TRUE),
              "row_id", "loc", "year")
    covariates <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], drop)
  }
  X <- as.matrix(d[, covariates, drop = FALSE])
  Y <- as.matrix(d[, resp, drop = FALSE])

  cells   <- list("00"=c(0,0), "10"=c(1,0), "01"=c(0,1), "11"=c(1,1))
  cellkey <- paste0(Z[,1], Z[,2])
  ck_lv   <- names(cells)

  ## plug-in outcome surfaces ghat_{S'} (n x q) from the posterior mean
  ghat <- lapply(cells, function(z) regime_outcome(fit, z, test = FALSE)$mean)
  names(ghat) <- ck_lv

  ## cross-fitted cell propensities ehat (n x 4)
  fit_prop <- function(Xtr, ytr, Xte) {
    if (propensity %in% c("auto","ranger") && requireNamespace("ranger", quietly = TRUE)) {
      rf <- ranger::ranger(y ~ ., data = data.frame(y = ytr, Xtr),
                           probability = TRUE, num.trees = 500)
      p <- stats::predict(rf, data.frame(Xte))$predictions
      return(p[, match(ck_lv, colnames(p)), drop = FALSE])
    }
    if (propensity %in% c("auto","multinom") && requireNamespace("nnet", quietly = TRUE)) {
      mn <- nnet::multinom(y ~ ., data = data.frame(y = ytr, Xtr), trace = FALSE, maxit = 300)
      p <- stats::predict(mn, newdata = data.frame(Xte), type = "probs")
      return(p[, match(ck_lv, colnames(p)), drop = FALSE])
    }
    sapply(ck_lv, function(cc) {                     # per-cell logistic fallback
      yb <- as.integer(ytr == cc)
      g  <- stats::glm(yb ~ ., data = data.frame(yb = yb, Xtr), family = stats::binomial())
      stats::predict(g, newdata = data.frame(Xte), type = "response")
    })
  }
  set.seed(seed); fold <- sample(rep(seq_len(folds), length.out = n))
  ehat <- matrix(NA_real_, n, length(cells), dimnames = list(NULL, ck_lv))
  for (f in seq_len(folds)) {
    te <- which(fold == f); trn <- which(fold != f)
    ehat[te, ] <- fit_prop(X[trn,,drop=FALSE],
                           factor(cellkey[trn], levels = ck_lv), X[te,,drop=FALSE])
  }
  ehat[is.na(ehat)] <- 1/length(cells)
  ehat <- pmin(pmax(ehat, trim), 1 - trim)

  ## one-step cell means chi (4 x q) and centered influence contributions phi (n x 4 x q)
  chi <- matrix(0, length(cells), q, dimnames = list(ck_lv, resp))
  phi <- array(0, c(n, length(cells), q), dimnames = list(NULL, ck_lv, resp))
  for (ci in seq_along(cells)) {
    cc <- ck_lv[ci]; m <- as.integer(cellkey == cc); g <- ghat[[cc]]; e <- ehat[, ci]
    for (k in seq_len(q)) {
      psi_i <- g[, k] + m / e * (Y[, k] - g[, k])
      chi[ci, k] <- mean(psi_i)
      phi[, ci, k] <- psi_i - chi[ci, k]
    }
  }

  ## Mobius combinations for the 2x2 estimands
  combos <- list(
    tau_Z1      = c("10"=1, "00"=-1),
    tau_Z2      = c("01"=1, "00"=-1),
    `tau_Z1:Z2` = c("11"=1, "10"=-1, "01"=-1, "00"=1),
    joint       = c("11"=1, "00"=-1))
  zq <- stats::qnorm(1 - (1 - level)/2)

  ## raw posterior plug-in ATE for comparison (component posterior means)
  raw_of <- function(nm, k) {
    if (nm == "joint")
      colMeans(contrast_effect(fit, rep(1, 2), rep(0, 2), test = FALSE)$mean)[k]
    else colMeans(component_effect(fit, nm, test = FALSE)$mean)[k]
  }

  rows <- list()
  for (nm in names(combos)) {
    w <- combos[[nm]]; ck <- names(w)
    for (k in seq_len(q)) {
      est  <- sum(w * chi[ck, k])
      sub  <- matrix(phi[, ck, k, drop = FALSE], nrow = n)  # n x |ck| (drop the k dim)
      infl <- as.numeric(sub %*% as.numeric(w))             # n-vector
      se   <- sqrt(stats::var(infl) / n)
      rows[[length(rows)+1]] <- data.frame(
        component = nm, outcome = resp[k],
        raw_plugin = round(raw_of(nm, k), 4),
        debiased = round(est, 4), se = round(se, 4),
        lo = round(est - zq*se, 4), hi = round(est + zq*se, 4),
        excludes_zero = (est - zq*se > 0) | (est + zq*se < 0),
        row.names = NULL)
    }
  }
  out <- do.call(rbind, rows)
  class(out) <- c("mvbcf_debiased", "data.frame")
  out
}

#' @export
print.mvbcf_debiased <- function(x, ...) {
  cat("Cross-fitted debiased average effects (one-step / AIPW; efficient SE, Wald CI)\n\n")
  print.data.frame(x, row.names = FALSE); invisible(x)
}
