#' Split records for one component's forest (var, value, iter, tree, comp).
#' @keywords internal
.splits_multi <- function(fit, component = "mu") {
  k <- if (is.numeric(component)) component else match(component, fit$components$labels)
  prog <- fit$components$is_prognostic[k]
  varNames <- if (prog) colnames(fit$X_con) else colnames(fit$X_mod)
  tdf <- fit$model$tree_df[[k]]
  sp <- as.data.frame(t(tdf))[-1, ]
  names(sp) <- c("var", "value", "iteration", "treeNum", "component")
  sp <- sp[sp$var != 0, ]
  sp$name <- varNames[sp$var]
  sp
}

#' Variable importance (inclusion proportion) for a component's forest
#'
#' Chipman et al. (2010) inclusion-proportion importance computed directly from
#' the stored split records of the chosen component (\code{"mu"}, \code{"tau_Z1"},
#' \code{"tau_Z1:Z2"}, ...). Robust to the many sparse genotype dummies that make
#' native bartMan importance unstable.
#'
#' @param fit A \code{\link{fit_mvbcf_multi}} object.
#' @param component Component label or index (default \code{"mu"}).
#' @param aggregate_gen Collapse the genotype dummies into one entry.
#' @return A data frame ordered by mean inclusion proportion.
#' @export
vimp_table <- function(fit, component = "mu", aggregate_gen = TRUE) {
  sp <- .splits_multi(fit, component)
  if (nrow(sp) == 0) return(data.frame())
  sp$grp <- if (aggregate_gen) ifelse(grepl("^gen_", sp$name), "genotype (aggregated)", sp$name) else sp$name
  tab <- table(sp$iteration, sp$grp); z <- sweep(tab, 1, rowSums(tab), "/")
  zdf <- as.data.frame.table(z, responseName = "z")
  out <- do.call(rbind, lapply(split(zdf, zdf$Var2), function(s)
    data.frame(variable = as.character(s$Var2[1]), propMean = mean(s$z),
               lowerQ = stats::quantile(s$z, .25), median = stats::median(s$z),
               upperQ = stats::quantile(s$z, .75), stringsAsFactors = FALSE)))
  out[order(-out$propMean), ]
}

#' Variable interactions (co-occurrence) within a component's trees
#' @inheritParams vimp_table
#' @param vars Variables to consider (default: the covariates).
#' @param top Keep the top-N pairs (\code{NULL} = all).
#' @export
vint_table <- function(fit, component = "mu", vars = fit$covariates, top = 10) {
  sp <- .splits_multi(fit, component); sp <- sp[sp$name %in% vars, ]
  if (nrow(sp) == 0) return(data.frame())
  key <- paste(sp$iteration, sp$treeNum, sep = "_")
  it <- tapply(sp$iteration, key, function(x) x[1])
  prs <- tapply(sp$name, key, function(nm) {
    u <- sort(unique(nm)); if (length(u) < 2) return(character(0))
    apply(utils::combn(u, 2), 2, paste, collapse = " : ")
  })
  df <- do.call(rbind, lapply(names(prs), function(k) {
    if (length(prs[[k]]) == 0) return(NULL)
    data.frame(iteration = it[[k]], pair = prs[[k]], stringsAsFactors = FALSE)
  }))
  if (is.null(df)) return(data.frame())
  tab <- table(df$iteration, df$pair); z <- sweep(tab, 1, rowSums(tab), "/")
  zdf <- as.data.frame.table(z, responseName = "z")
  out <- do.call(rbind, lapply(split(zdf, zdf$Var2), function(s)
    data.frame(pair = as.character(s$Var2[1]), propMean = mean(s$z),
               lowerQ = stats::quantile(s$z, .25), upperQ = stats::quantile(s$z, .75),
               stringsAsFactors = FALSE)))
  out <- out[order(-out$propMean), ]
  if (!is.null(top)) out <- utils::head(out, top)
  out
}

#' Plot inclusion-proportion importance for a component
#' @inheritParams vimp_table
#' @export
plot_vimp <- function(fit, component = "mu", aggregate_gen = TRUE) {
  v <- vimp_table(fit, component, aggregate_gen)
  v$variable <- factor(v$variable, levels = rev(v$variable))
  ggplot2::ggplot(v, ggplot2::aes(propMean, variable)) +
    ggplot2::geom_segment(ggplot2::aes(x = lowerQ, xend = upperQ, yend = variable), colour = "grey60") +
    ggplot2::geom_point(size = 2, colour = "#2166AC") +
    ggplot2::labs(title = paste0("Variable importance - ", if(is.numeric(component)) fit$components$labels[component] else component),
                  subtitle = "point = mean; bar = 25-75% quantile interval",
                  x = "inclusion proportion", y = NULL) +
    ggplot2::theme_minimal(base_size = 10)
}

#' Build a bartMan trees_data object from a component's forest
#'
#' Lets you use native \pkg{bartMan} plotting (e.g. \code{plotTrees},
#' \code{splitDensity}) on any component forest of the multi-treatment MVBCF.
#'
#' @param fit A fit object.
#' @param component Component label/index (default "mu").
#' @return A bartMan trees data object.
#' @export
mvbcf_trees <- function(fit, component = "mu") {
  if (!requireNamespace("bartMan", quietly = TRUE)) stop("please install 'bartMan'")
  k <- if (is.numeric(component)) component else match(component, fit$components$labels)
  prog <- fit$components$is_prognostic[k]
  Xdes <- if (prog) fit$X_con else fit$X_mod
  f_data <- data.frame(y = fit$data[[fit$responses[1]]], Xdes, check.names = FALSE)
  tr <- t(fit$model$tree_df[[k]]); tr <- tr[-1, 1:4]
  tr[, 1] <- ifelse(tr[, 1] == 0, NA, tr[, 1])
  colnames(tr) <- c("var", "value", "iteration", "treeNum")
  bartMan::tree_dataframe(data = f_data, trees = tr, response = "y")
}
