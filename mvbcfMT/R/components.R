#' Build the additive-component specification for a multi-treatment MVBCF.
#'
#' Given the treatment columns and a maximum interaction order \code{r}, this
#' enumerates the ANOVA/Mobius components of the factorial response surface:
#' the prognostic term \eqn{\mu} (the empty subset) plus one forest per subset
#' \eqn{S \subseteq \{1,\dots,T\}} with \eqn{1 \le |S| \le r}. For each component
#' it returns the indicator vector \eqn{D^{(S)}_i = \prod_{t\in S} Z_{it}} and a
#' human-readable label.
#'
#' @param Zmat An \code{n x T} 0/1 matrix of treatment indicators.
#' @param treat_names Length-T character vector of treatment names.
#' @param order Maximum interaction order to include (1 = main effects only,
#'   2 = up to pairwise, ... , T = full factorial).
#' @return A list with \code{labels}, \code{subsets} (list of integer index
#'   vectors; \code{integer(0)} = prognostic), \code{indicator} (n x K matrix)
#'   and \code{is_prognostic} (logical K).
#' @keywords internal
#' @export
build_components <- function(Zmat, treat_names, order = 2L) {
  Ti <- ncol(Zmat)
  order <- min(order, Ti)
  subsets <- list(integer(0))                     # prognostic = empty subset
  for (o in seq_len(order))
    subsets <- c(subsets, utils::combn(Ti, o, simplify = FALSE))
  labels <- vapply(subsets, function(S) {
    if (length(S) == 0) "mu"
    else paste0("tau_", paste(treat_names[S], collapse = ":"))
  }, character(1))
  ind <- vapply(subsets, function(S) {
    if (length(S) == 0) rep(1, nrow(Zmat))
    else apply(Zmat[, S, drop = FALSE], 1, prod)
  }, numeric(nrow(Zmat)))
  ind <- matrix(ind, nrow = nrow(Zmat))
  list(labels = labels, subsets = subsets, indicator = ind,
       is_prognostic = vapply(subsets, function(S) length(S) == 0, logical(1)))
}

#' Replicate a length-n indicator across q outcome columns (engine expects n x q).
#' @keywords internal
.rep_q <- function(z, q) matrix(rep(z, q), ncol = q)

#' Estimate a propensity score for a binary treatment.
#' Uses \pkg{dbarts} if available (matching the original MVBCF), otherwise a
#' logistic GLM on the design (adequate and dependency-free for demonstration).
#' @keywords internal
.propensity <- function(X, z, seed = 1) {
  if (requireNamespace("dbarts", quietly = TRUE)) {
    set.seed(seed)
    pm <- dbarts::bart(x.train = X, y.train = z, verbose = FALSE)
    return(colMeans(stats::pnorm(pm$yhat.train)))
  }
  df <- data.frame(z = z, X)
  stats::fitted(stats::glm(z ~ ., data = df, family = stats::binomial()))
}
