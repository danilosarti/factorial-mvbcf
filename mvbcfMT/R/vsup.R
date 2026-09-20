# Interpretability layer: value-suppressing uncertainty palette (VSUP) and
# per-component importance, computed directly from the stored posterior (no
# dependency on bartMan/vivid, which are optional and unstable on the sparse
# genotype dummies -- same robust route as the mvbcfMET pipeline).

#' Value-Suppressing Uncertainty Palette: map (value, uncertainty) -> hex.
#' Hue/lightness encodes the value on a diverging ramp; rising uncertainty
#' desaturates toward grey, so uncertain cells "recede" (Correll et al. 2018).
#' @param value numeric vector (e.g. posterior-mean effect).
#' @param uncertainty non-negative vector (e.g. posterior SD or CI half-width).
#' @param palette diverging colours for low/mid/high value.
#' @param mid central value mapped to the neutral colour (default 0).
#' @param max_desat maximum fraction blended to grey at the highest uncertainty.
#' @param grey the grey blended toward.
#' @return character vector of hex colours.
#' @export
#' @keywords internal
vsup_palette <- function(value, uncertainty,
                         palette = c("#B2182B", "#F7F7F7", "#2166AC"),
                         mid = 0, max_desat = 0.85, grey = "#BDBDBD") {
  v <- value - mid
  rng <- max(abs(v), na.rm = TRUE); if (rng == 0 || !is.finite(rng)) rng <- 1
  t <- pmin(pmax((v / rng + 1) / 2, 0), 1)                 # -> [0,1], 0.5 = mid
  ramp <- grDevices::colorRamp(palette, space = "Lab")
  base <- ramp(t) / 255                                     # n x 3 in [0,1]
  u <- uncertainty; urng <- stats::quantile(u, 0.95, na.rm = TRUE)
  if (!is.finite(urng) || urng == 0) urng <- max(u, na.rm = TRUE); if (urng == 0) urng <- 1
  w <- pmin(u / urng, 1) * max_desat                        # blend weight to grey
  g <- grDevices::col2rgb(grey)[, 1] / 255
  mixed <- base * (1 - w) + matrix(g, nrow(base), 3, byrow = TRUE) * w
  grDevices::rgb(mixed[, 1], mixed[, 2], mixed[, 3])
}

#' G x E tile of a treatment effect with a value-suppressing uncertainty palette
#'
#' Renders one component's per-cell effect on a genotype x environment tile,
#' colouring by the posterior-mean effect but \emph{suppressing} (greying) cells
#' whose posterior is uncertain. Cells never grown in the field (\code{observed==0})
#' are marked. Requires a fit built with \code{gen}/\code{env}.
#'
#' @param fit A \code{\link{fit_mvbcf_multi}} object with a G x E grid.
#' @param component Component label (e.g. \code{"tau_Z1"}, \code{"tau_Z1:Z2"}).
#' @param response Which response column to show.
#' @param max_desat,mid See \code{\link{vsup_palette}}.
#' @return A ggplot.
#' @export
plot_effect_vsup <- function(fit, component, response = fit$responses[1],
                             max_desat = 0.85, mid = 0) {
  if (is.null(fit$grid)) stop("this fit has no G x E grid (build with gen=, env=)")
  k <- if (is.numeric(component)) component else match(component, fit$components$labels)
  cube <- fit$model$predictions_test[[k]]
  ri <- match(response, fit$responses)
  m  <- apply(cube[, ri, , drop = TRUE], 1, mean)
  sdv <- apply(cube[, ri, , drop = TRUE], 1, stats::sd)
  df <- data.frame(fit$grid[, c(fit$gen, fit$env, "observed")],
                   eff = m, unc = sdv)
  df$col <- vsup_palette(df$eff, df$unc, mid = mid, max_desat = max_desat)
  names(df)[1:2] <- c("gen", "env")
  ggplot2::ggplot(df, ggplot2::aes(env, gen)) +
    ggplot2::geom_tile(ggplot2::aes(fill = col), colour = "white", linewidth = 0.2) +
    ggplot2::geom_point(data = df[df$observed == 0, ],
                        shape = 4, size = 0.7, colour = "grey30", alpha = 0.5) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(title = paste0("VSUP: ", fit$components$labels[k], "  (", response, ")"),
                  subtitle = "hue = effect; greyed = uncertain; x = never grown",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 5),
                   panel.grid = ggplot2::element_blank())
}

#' A standalone VSUP legend (value x uncertainty grid)
#' @param mid,max_desat See \code{\link{vsup_palette}}.
#' @param value_lab,unc_lab axis labels.
#' @return A ggplot of the 2-D legend.
#' @export
vsup_legend <- function(mid = 0, max_desat = 0.85,
                        value_lab = "effect", unc_lab = "uncertainty (posterior SD)") {
  g <- expand.grid(value = seq(-1, 1, length.out = 40),
                   unc = seq(0, 1, length.out = 20))
  g$col <- vsup_palette(g$value, g$unc, mid = mid, max_desat = max_desat)
  ggplot2::ggplot(g, ggplot2::aes(value, unc, fill = col)) +
    ggplot2::geom_raster() + ggplot2::scale_fill_identity() +
    ggplot2::labs(title = "VSUP legend", x = value_lab, y = unc_lab) +
    ggplot2::theme_minimal(base_size = 9)
}

#' Context-agnostic VSUP over any two covariates
#'
#' Domain-neutral alternative to \code{\link{plot_effect_vsup}} for fits without a
#' genotype x environment grid: bins two chosen covariates into a 2-D grid, and
#' tiles the average effect of a component coloured by the value with uncertainty
#' suppression (greying where the posterior is diffuse or the cell is sparse).
#'
#' @param fit A fit object.
#' @param component Component label/index.
#' @param xvar,yvar Covariate column names to place on the axes.
#' @param response Which response.
#' @param bins Number of bins per axis.
#' @param test Use test rows (default FALSE).
#' @return A ggplot.
#' @export
plot_effect_vsup_xy <- function(fit, component, xvar, yvar,
                                response = fit$responses[1], bins = 8, test = FALSE) {
  k <- if (is.numeric(component)) component else match(component, fit$components$labels)
  cube <- if (test) fit$model$predictions_test[[k]] else fit$model$predictions[[k]]
  ri <- match(response, fit$responses)
  m  <- apply(cube[, ri, , drop = TRUE], 1, mean)
  sdv <- apply(cube[, ri, , drop = TRUE], 1, stats::sd)
  X <- if (test && !is.null(fit$grid)) fit$grid else fit$data
  xb <- cut(X[[xvar]], breaks = bins); yb <- cut(X[[yvar]], breaks = bins)
  agg <- stats::aggregate(data.frame(eff = m, unc = sdv, n = 1),
                   by = list(x = xb, y = yb), FUN = function(z) z)
  ag <- do.call(data.frame, lapply(list(
    x = agg$x, y = agg$y,
    eff = tapply(m, list(xb, yb), mean)[cbind(agg$x, agg$y)],
    unc = tapply(sdv, list(xb, yb), mean)[cbind(agg$x, agg$y)],
    n   = as.integer(table(xb, yb)[cbind(agg$x, agg$y)])), identity))
  ag <- ag[stats::complete.cases(ag), ]
  # inflate uncertainty where cells are sparse
  ag$unc2 <- ag$unc / pmax(sqrt(ag$n), 1) * max(sqrt(ag$n))
  ag$col <- vsup_palette(ag$eff, ag$unc2, mid = 0)
  ggplot2::ggplot(ag, ggplot2::aes(x, y, fill = col)) +
    ggplot2::geom_tile(colour = "white") + ggplot2::scale_fill_identity() +
    ggplot2::labs(title = paste0("VSUP: ", fit$components$labels[k], " (", response, ")"),
                  subtitle = "hue = effect; greyed = uncertain / sparse cell",
                  x = xvar, y = yvar) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank())
}

#' Forest plot of the average component effects (with 95\% CrI)
#'
#' Context-agnostic summary: the posterior mean and 95\% credible interval of the
#' average effect of every component, for one response. Works in any domain.
#'
#' @param fit A fit object.
#' @param response Which response.
#' @param test Use test rows/grid (default FALSE = training rows).
#' @return A ggplot.
#' @export
plot_effect_forest <- function(fit, response = fit$responses[1], test = FALSE) {
  ri <- match(response, fit$responses)
  labs <- fit$components$labels
  labs <- labs[!fit$components$is_prognostic]                 # effects only
  dat <- do.call(rbind, lapply(labs, function(l) {
    ce <- component_effect(fit, l, test = test)
    data.frame(component = l, ate = mean(ce$mean[, ri]),
               lo = mean(ce$lo[, ri]), hi = mean(ce$hi[, ri]))
  }))
  dat$component <- factor(dat$component, levels = rev(dat$component))
  ggplot2::ggplot(dat, ggplot2::aes(ate, component)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = .2, colour = "grey50") +
    ggplot2::geom_point(colour = "#B2182B", size = 2.4) +
    ggplot2::labs(title = paste0("Average effects on ", response, " (95% CrI)"),
                  x = "effect", y = NULL) +
    ggplot2::theme_minimal(base_size = 10)
}

#' Sequential value-suppressing uncertainty palette (for non-negative values)
#' @param value non-negative vector.
#' @param uncertainty non-negative vector.
#' @param low,high sequential endpoints.
#' @param max_desat,grey suppression.
#' @param power Exponent applied to the normalised value ramp (default 0.6).
#' @export
vsup_palette_seq <- function(value, uncertainty, low = "#F7FBFF", high = "#08306B",
                             max_desat = 0.75, grey = "#C7C7C7", power = 0.6) {
  v <- value; rng <- max(v, na.rm = TRUE); if (!is.finite(rng) || rng == 0) rng <- 1
  t <- pmin(pmax(v / rng, 0), 1)^power
  base <- grDevices::colorRamp(c(low, high), space = "Lab")(t) / 255
  u <- uncertainty; urng <- stats::quantile(u, 0.95, na.rm = TRUE)
  if (!is.finite(urng) || urng == 0) urng <- max(u, na.rm = TRUE); if (urng == 0) urng <- 1
  w <- pmin(u / urng, 1) * max_desat
  g <- grDevices::col2rgb(grey)[, 1] / 255
  mixed <- base * (1 - w) + matrix(g, nrow(base), 3, byrow = TRUE) * w
  grDevices::rgb(mixed[, 1], mixed[, 2], mixed[, 3])
}

#' VIVI matrix (variable importance + interaction) with posterior uncertainty
#'
#' For one component's forest, returns the VIVI matrix in the sense of the
#' \pkg{vivid}/\pkg{bartMan} value-suppressing displays: the diagonal holds each
#' variable's inclusion-proportion importance and the off-diagonal holds pairwise
#' within-tree interaction proportions (Chipman et al.\ 2010). Crucially it also
#' returns the \emph{posterior uncertainty} of every entry --- the standard
#' deviation across the stored MCMC iterations --- computed directly from the
#' per-iteration split records (no \pkg{bartMan} dependency). Fit with a longer
#' \code{tree_iters} for a well-estimated uncertainty.
#'
#' @param fit A \code{\link{fit_mvbcf_multi}} object.
#' @param component Component label/index (e.g. \code{"tau_Z1"}).
#' @param vars Variables to include (default: the covariates).
#' @return A list with \code{value} and \code{uncertainty} matrices (vars x vars).
#' @export
vivi_matrix <- function(fit, component = "mu", vars = fit$covariates) {
  sp <- .splits_multi(fit, component); sp <- sp[sp$name %in% vars, ]
  V <- length(vars); idx <- stats::setNames(seq_len(V), vars)
  iters <- sort(unique(sp$iteration))
  imp <- matrix(0, length(iters), V)
  intr <- array(0, dim = c(length(iters), V, V))
  for (i in seq_along(iters)) {
    si <- sp[sp$iteration == iters[i], ]
    if (nrow(si) == 0) next
    tab <- table(factor(si$name, levels = vars))
    imp[i, ] <- as.numeric(tab) / sum(tab)
    trees <- unique(si$treeNum); ntree <- length(trees)
    for (tn in trees) {
      nm <- unique(si$name[si$treeNum == tn]); if (length(nm) < 2) next
      cb <- utils::combn(nm, 2)
      for (c in seq_len(ncol(cb))) {
        a <- idx[cb[1, c]]; b <- idx[cb[2, c]]
        intr[i, a, b] <- intr[i, a, b] + 1; intr[i, b, a] <- intr[i, b, a] + 1
      }
    }
    if (ntree > 0) intr[i, , ] <- intr[i, , ] / ntree
  }
  Mval <- Munc <- matrix(0, V, V, dimnames = list(vars, vars))
  diag(Mval) <- colMeans(imp); diag(Munc) <- apply(imp, 2, stats::sd)
  for (a in seq_len(V)) for (b in seq_len(V)) if (a != b) {
    Mval[a, b] <- mean(intr[, a, b]); Munc[a, b] <- stats::sd(intr[, a, b])
  }
  list(value = Mval, uncertainty = Munc)
}

#' VSUP heatmap of a component's VIVI matrix (value suppressed by uncertainty)
#'
#' Renders \code{\link{vivi_matrix}} as a heatmap where hue encodes importance /
#' interaction and cells with uncertain posteriors recede to grey --- the
#' value-suppressing uncertainty display of \pkg{vivid}/\pkg{bartMan}, computed here
#' directly from the stored posterior so it applies to \emph{each} effect forest
#' (\eqn{\tau_1,\tau_2,\tau_{12}}), not just a single BART model.
#'
#' @param fit A fit object.
#' @param component Component label/index.
#' @param vars Variables (default covariates).
#' @param title Optional title.
#' @return A ggplot.
#' @export
plot_vivi_vsup <- function(fit, component = "mu", vars = fit$covariates, title = NULL) {
  vm <- vivi_matrix(fit, component, vars)
  lab <- if (is.numeric(component)) fit$components$labels[component] else component
  df <- expand.grid(row = vars, col = vars, stringsAsFactors = FALSE)
  df$value <- mapply(function(r, c) vm$value[r, c], df$row, df$col)
  df$unc   <- mapply(function(r, c) vm$uncertainty[r, c], df$row, df$col)
  df$col_hex <- vsup_palette_seq(df$value, df$unc)
  df$row <- factor(df$row, levels = rev(vars)); df$col <- factor(df$col, levels = vars)
  ggplot2::ggplot(df, ggplot2::aes(col, row, fill = col_hex)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(title = if (is.null(title)) paste0("VIVI (VSUP): ", lab) else title,
                  subtitle = "diag = importance, off-diag = interaction; greyed = uncertain",
                  x = NULL, y = NULL) +
    ggplot2::coord_equal() + ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank())
}

#' 2-D sequential VSUP legend (value x uncertainty)
#' @export
vsup_legend_seq <- function() {
  g <- expand.grid(value = seq(0, 1, length.out = 40), unc = seq(0, 1, length.out = 20))
  g$col <- vsup_palette_seq(g$value, g$unc)
  ggplot2::ggplot(g, ggplot2::aes(value, unc, fill = col)) +
    ggplot2::geom_raster() + ggplot2::scale_fill_identity() +
    ggplot2::labs(title = "VSUP legend", x = "importance / interaction",
                  y = "posterior uncertainty") +
    ggplot2::theme_minimal(base_size = 8)
}

#' @keywords internal
.vsup_mix <- function(t, w, low, high, grey, power = 0.7) {
  t <- pmin(pmax(t, 0), 1)^power
  base <- grDevices::colorRamp(c(low, high), space = "Lab")(t) / 255
  g <- grDevices::col2rgb(grey)[, 1] / 255
  m <- base * (1 - w) + matrix(g, nrow(base), 3, byrow = TRUE) * w
  grDevices::rgb(m[, 1], m[, 2], m[, 3])
}

#' vivid-style VIVI heatmap with value-suppressing uncertainty (two colour scales)
#'
#' Reproduces the \pkg{vivid}/\pkg{bartMan} value-suppressing display without those
#' packages: a single heatmap of one component's forest in which the diagonal
#' encodes variable \emph{importance} on a blue scale (Vimp) and the off-diagonal
#' encodes pairwise \emph{interaction} on a red scale (Vint), each suppressed toward
#' a warm grey as the posterior coefficient of variation (CV = SD/value across MCMC
#' iterations) grows. Pair with \code{\link{vsup_wedge}} for the fan legends.
#'
#' @param fit A fit object.
#' @param component Component label/index.
#' @param vars Variables on the axes.
#' @param cvmax CV at which cells are fully greyed.
#' @param max_desat Maximum grey blend.
#' @param title Optional title.
#' @return A ggplot; its attribute \code{"maxes"} holds the importance/interaction
#'   maxima for the matching legends.
#' @export
plot_vivi_vivid <- function(fit, component = "mu", vars = fit$covariates,
                            cvmax = 2, max_desat = 0.85, title = NULL) {
  vm <- vivi_matrix(fit, component, vars); V <- length(vars)
  lab <- if (is.numeric(component)) fit$components$labels[component] else component
  df <- expand.grid(ri = seq_len(V), ci = seq_len(V))
  df$value <- mapply(function(r, c) vm$value[r, c], df$ri, df$ci)
  df$sd    <- mapply(function(r, c) vm$uncertainty[r, c], df$ri, df$ci)
  df$diag  <- df$ri == df$ci
  df$cv <- pmin(df$sd / pmax(df$value, 1e-6), cvmax) / cvmax
  df$w  <- df$cv * max_desat
  vmax_d <- max(df$value[df$diag]); vmax_o <- max(df$value[!df$diag], na.rm = TRUE)
  di <- df$diag; oi <- !df$diag; df$hex <- NA_character_
  df$hex[di] <- .vsup_mix(df$value[di] / max(vmax_d, 1e-9), df$w[di], "#FFFFCC", "#08306B", "#BFB8AE")
  df$hex[oi] <- .vsup_mix(df$value[oi] / max(vmax_o, 1e-9), df$w[oi], "#FFFFCC", "#67000D", "#BFB8AE")
  df$row <- factor(vars[df$ri], levels = rev(vars)); df$col <- factor(vars[df$ci], levels = vars)
  g <- ggplot2::ggplot(df, ggplot2::aes(col, row, fill = hex)) +
    ggplot2::geom_tile(colour = "grey88", linewidth = 0.3) +
    ggplot2::scale_fill_identity() + ggplot2::coord_equal() +
    ggplot2::labs(title = if (is.null(title)) paste0("VIVI (VSUP): ", lab) else title,
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 6), panel.grid = ggplot2::element_blank())
  attr(g, "maxes") <- c(imp = vmax_d, int = vmax_o); g
}

#' Fan (wedge) VSUP legend, in the \pkg{vivid} style
#'
#' Draws a quarter-fan legend: angular position encodes the value (importance or
#' interaction) and radial distance encodes the posterior CV; colour is the same
#' value-suppressing map as \code{\link{plot_vivi_vivid}}.
#'
#' @param vmax Maximum value on the value (angular) axis.
#' @param kind \code{"imp"} (blue) or \code{"int"} (red).
#' @param cvmax Maximum CV on the radial axis.
#' @param nv,nu Value/CV bands.
#' @param title Legend title.
#' @param max_desat Maximum fraction blended toward grey at the highest CV.
#' @return A ggplot.
#' @export
vsup_wedge <- function(vmax, kind = c("int", "imp"), cvmax = 2, nv = 6, nu = 5,
                       title = NULL, max_desat = 0.85) {
  kind <- match.arg(kind); high <- if (kind == "int") "#67000D" else "#08306B"
  phi <- seq(pi * 0.78, pi * 0.22, length.out = nv + 1)   # value 0 -> max (left->right)
  rho <- seq(0.28, 1.0, length.out = nu + 1)              # CV 0 -> cvmax (in->out)
  polys <- list(); k <- 0
  for (i in seq_len(nv)) for (j in seq_len(nu)) { k <- k + 1
    p1 <- phi[i]; p2 <- phi[i + 1]; r1 <- rho[j]; r2 <- rho[j + 1]
    xs <- c(r1 * cos(p1), r1 * cos(p2), r2 * cos(p2), r2 * cos(p1))
    ys <- c(r1 * sin(p1), r1 * sin(p2), r2 * sin(p2), r2 * sin(p1))
    val <- (i - 0.5) / nv; w <- ((j - 0.5) / nu) * max_desat
    polys[[k]] <- data.frame(id = k, x = xs, y = ys,
                             hex = .vsup_mix(val, w, "#FFFFCC", high, "#BFB8AE"))
  }
  d <- do.call(rbind, polys)
  # tick labels: value along outer arc, CV along right radius
  vt <- data.frame(v = seq(0, vmax, length.out = nv + 1),
                   x = 1.14 * cos(phi), y = 1.14 * sin(phi))
  ct <- data.frame(c = round(seq(0, cvmax, length.out = nu + 1), 2),
                   x = rho * cos(phi[nv + 1]) + 0.06, y = rho * sin(phi[nv + 1]))
  ggplot2::ggplot(d, ggplot2::aes(x, y, group = id, fill = hex)) +
    ggplot2::geom_polygon(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_identity() + ggplot2::coord_equal() +
    ggplot2::geom_text(data = vt, ggplot2::aes(x, y, label = signif(v, 2)),
                       inherit.aes = FALSE, size = 2.3) +
    ggplot2::geom_text(data = ct, ggplot2::aes(x, y, label = c),
                       inherit.aes = FALSE, size = 2.1, hjust = 0) +
    ggplot2::labs(title = if (is.null(title)) (if (kind == "int") "Vint" else "Vimp") else title,
                  subtitle = "arc = value; radius = CV") +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Faceted per-component variable importance
#'
#' Inclusion-proportion importance for several components at once (\eqn{\mu} and
#' each effect forest), so one can read off which covariates drive the baseline
#' vs each treatment effect vs the interaction.
#'
#' @param fit A fit object.
#' @param components Which components (default: all).
#' @param top Keep the top-N variables per component.
#' @param aggregate_gen Collapse genotype dummies.
#' @return A ggplot.
#' @export
plot_component_importance <- function(fit, components = fit$components$labels,
                                      top = 8, aggregate_gen = TRUE) {
  dat <- do.call(rbind, lapply(components, function(cl) {
    v <- vimp_table(fit, cl, aggregate_gen = aggregate_gen)
    if (nrow(v) == 0) return(NULL)
    v <- utils::head(v, top); v$component <- cl; v
  }))
  dat$component <- factor(dat$component, levels = components)
  dat$variable <- stats::reorder(dat$variable, dat$propMean)
  ggplot2::ggplot(dat, ggplot2::aes(propMean, variable)) +
    ggplot2::geom_segment(ggplot2::aes(x = lowerQ, xend = upperQ, yend = variable), colour = "grey65") +
    ggplot2::geom_point(colour = "#2166AC", size = 1.6) +
    ggplot2::facet_wrap(~component, scales = "free", ncol = 2) +
    ggplot2::labs(title = "Variable importance by component (inclusion proportion)",
                  x = "inclusion proportion", y = NULL) +
    ggplot2::theme_minimal(base_size = 9)
}
