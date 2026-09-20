#' Fit a multi-treatment multivariate Bayesian Causal Forest
#'
#' Fits the additive factorial MVBCF
#' \deqn{Y_i = \mu(x_i) + \sum_{\emptyset \ne S,\,|S|\le r}
#'        \Big(\prod_{t\in S} Z_{it}\Big)\, \tau_S(x_i) + \varepsilon_i,
#'        \quad \varepsilon_i \sim N_q(0,\Sigma),}
#' where each \eqn{\mu} and \eqn{\tau_S} is its own multivariate sum-of-trees
#' forest. With two treatments and \code{order = 2} this is the
#' \eqn{(\mu,\tau_1,\tau_2,\tau_{12})} factorial model; \code{order} truncates
#' the interaction lattice for more than two treatments.
#'
#' @param data A data frame, one row per unit.
#' @param responses Character vector of response columns (q >= 1).
#' @param treatments Character vector of binary (0/1) treatment columns.
#' @param covariates Moderator/control columns. If \code{NULL}, all numeric
#'   columns except responses, treatments, ids and any \code{true_*}/ground-truth
#'   columns are used.
#' @param order Maximum interaction order (default 2 = up to pairwise).
#' @param gen,env Optional genotype / environment id columns (MET use). When
#'   given, genotype is one-hot encoded (prefix \code{gen_}) into the moderators
#'   and a full gen x env prediction grid is built.
#' @param n_iter,n_burn,keep_every,n_tree,n_tree_tau,min_nodesize MCMC/forest
#'   settings. \code{n_tree} is the prognostic forest size; \code{n_tree_tau} is
#'   used for every effect forest.
#' @param sd_control,sd_effect Prior scale multipliers: node prior covariance is
#'   \code{diag(sd_control^2 / n_tree)} for \eqn{\mu} and
#'   \code{diag(sd_effect^2 / n_tree_tau)} for every effect forest (Hahn-style
#'   half-scale on effects).
#' @param tree_iters Iterations whose tree structure is stored (bartMan). Default
#'   last 20.
#' @param seed Seed for the propensity models.
#' @param verbose Print C++ progress.
#' @return An object of class \code{"mvbcfMT"}.
#' @examples
#' \donttest{
#' sim <- simulate_multi(n = 150, q = 2, Ti = 2, seed = 1)
#' fit <- fit_mvbcf_multi(sim$data, responses = sim$responses,
#'                        treatments = sim$treatments, covariates = sim$covariates,
#'                        n_iter = 100, n_burn = 50, n_tree = 30, n_tree_tau = 20)
#' print(fit)
#' }
#' @export
fit_mvbcf_multi <- function(data, responses, treatments,
                            covariates = NULL, order = 2L,
                            gen = NULL, env = NULL,
                            n_iter = 1500, n_burn = 750, keep_every = 2,
                            n_tree = 100, n_tree_tau = 50, min_nodesize = 5,
                            sd_control = 1, sd_effect = 0.5,
                            tree_iters = NULL, seed = 1, verbose = FALSE) {
  stopifnot(length(treatments) >= 1)
  if (is.null(covariates)) {
    drop <- c(responses, treatments, gen, env, "row_id", "loc", "year",
              grep("^true_", names(data), value = TRUE), "severity", "tolerance")
    num <- names(data)[vapply(data, is.numeric, logical(1))]
    covariates <- setdiff(num, drop)
  }
  q <- length(responses)
  y <- as.matrix(data[, responses, drop = FALSE])
  Zmat <- as.matrix(data[, treatments, drop = FALSE]); storage.mode(Zmat) <- "double"

  # ---- moderator design (covariates, plus genotype one-hot for MET) ----------
  gen_levels <- env_levels <- NULL
  if (!is.null(gen)) {
    gen_levels <- sort(unique(data[[gen]]))
    oh <- sapply(gen_levels, function(l) as.integer(data[[gen]] == l))
    if (is.null(dim(oh))) oh <- matrix(oh, nrow = nrow(data))
    colnames(oh) <- paste0("gen_", gen_levels)
    X_mod <- cbind(oh, as.matrix(data[, covariates, drop = FALSE]))
  } else {
    X_mod <- as.matrix(data[, covariates, drop = FALSE])
  }
  if (!is.null(env)) env_levels <- sort(unique(data[[env]]))

  # ---- component specification (mu + effect forests) -------------------------
  comp <- build_components(Zmat, treatments, order = order)
  K <- length(comp$labels)

  # ---- propensity scores as prognostic controls ------------------------------
  ps <- vapply(seq_along(treatments), function(j)
    .propensity(X_mod, Zmat[, j], seed = seed + j), numeric(nrow(data)))
  colnames(ps) <- paste0("pro_", treatments)
  X_con <- cbind(X_mod, ps)

  # ---- test design: gen x env grid (MET) or the training design itself -------
  if (!is.null(gen) && !is.null(env)) {
    grid <- expand.grid(gen_levels, env_levels, stringsAsFactors = FALSE)
    names(grid) <- c(gen, env)
    env_tab <- unique(data[, c(env, covariates)])
    grid <- merge(grid, env_tab, by = env, all.x = TRUE)
    grid <- grid[, c(gen, env, covariates)]
    obs_key <- unique(paste(data[[gen]], data[[env]]))
    grid$observed <- as.integer(paste(grid[[gen]], grid[[env]]) %in% obs_key)
    oh_t <- sapply(gen_levels, function(l) as.integer(grid[[gen]] == l))
    if (is.null(dim(oh_t))) oh_t <- matrix(oh_t, nrow = nrow(grid))
    colnames(oh_t) <- paste0("gen_", gen_levels)
    X_mod_test <- cbind(oh_t, as.matrix(grid[, covariates, drop = FALSE]))
    ps_test <- matrix(rep(colMeans(ps), each = nrow(grid)), nrow = nrow(grid))
    colnames(ps_test) <- colnames(ps)
    X_con_test <- cbind(X_mod_test, ps_test)
  } else {
    grid <- NULL
    X_mod_test <- X_mod; X_con_test <- X_con
  }

  # ---- assemble engine inputs (K components) ---------------------------------
  X_list <- X_test_list <- D_list <- sp_list <- vector("list", K)
  a_vec <- b_vec <- nt_vec <- numeric(K); addm <- logical(K)
  for (k in seq_len(K)) {
    prog <- comp$is_prognostic[k]
    X_list[[k]]      <- if (prog) X_con      else X_mod
    X_test_list[[k]] <- if (prog) X_con_test else X_mod_test
    D_list[[k]]      <- .rep_q(comp$indicator[, k], q)
    nt_vec[k]        <- if (prog) n_tree else n_tree_tau
    sp_list[[k]]     <- if (prog) diag(sd_control^2 / n_tree, q)
                        else       diag(sd_effect^2  / n_tree_tau, q)
    a_vec[k]         <- if (prog) 0.95 else 0.25
    b_vec[k]         <- if (prog) 2    else 3
    addm[k]          <- prog
  }
  if (is.null(tree_iters)) tree_iters <- seq.int(max(1, n_iter - 19), n_iter)

  run <- function()
    fast_bart_multi(y = y, X_list = X_list, X_test_list = X_test_list,
                    D_list = D_list, alpha = a_vec, beta = b_vec,
                    sigma_par_list = sp_list, n_tree = nt_vec, add_mean = addm,
                    v_0 = 5, sigma_0 = diag(q),
                    n_iter = n_iter, min_nodesize = min_nodesize,
                    n_burn = n_burn, keep_every = keep_every,
                    tree_iters = as.numeric(tree_iters))
  if (verbose) model <- run()
  else invisible(utils::capture.output(model <- run()))

  structure(list(model = model, data = data, grid = grid, responses = responses,
                 treatments = treatments, covariates = covariates,
                 gen = gen, env = env, gen_levels = gen_levels, env_levels = env_levels,
                 components = comp, X_con = X_con, X_mod = X_mod,
                 propensity = ps, order = order,
                 settings = list(n_iter = n_iter, n_burn = n_burn, keep_every = keep_every,
                                 n_tree = n_tree, n_tree_tau = n_tree_tau,
                                 min_nodesize = min_nodesize, tree_iters = tree_iters)),
            class = "mvbcfMT")
}

#' Fit the two-treatment MVBCF (convenience wrapper)
#'
#' Thin wrapper around \code{\link{fit_mvbcf_multi}} for exactly two treatments,
#' fitting the full \eqn{(\mu,\tau_1,\tau_2,\tau_{12})} factorial model.
#'
#' @inheritParams fit_mvbcf_multi
#' @param treatments Length-2 character vector of the two treatment columns.
#' @param ... Additional arguments passed to \code{\link{fit_mvbcf_multi}}.
#' @export
fit_mvbcf2 <- function(data, responses, treatments, covariates = NULL, ...) {
  stopifnot(length(treatments) == 2)
  fit_mvbcf_multi(data, responses = responses, treatments = treatments,
                  covariates = covariates, order = 2L, ...)
}

#' @export
print.mvbcfMT <- function(x, ...) {
  cat("<mvbcfMT fit: multi-treatment multivariate BCF>\n")
  cat(" observations :", nrow(x$data), "\n")
  cat(" responses    :", paste(x$responses, collapse = ", "), "\n")
  cat(" treatments   :", paste(x$treatments, collapse = ", "),
      "  (interaction order <=", x$order, ")\n")
  cat(" components   :", paste(x$components$labels, collapse = ", "), "\n")
  if (!is.null(x$grid))
    cat(" G x E grid   :", nrow(x$grid), "cells (", sum(x$grid$observed), "observed )\n")
  invisible(x)
}
