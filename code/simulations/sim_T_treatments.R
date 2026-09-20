# =============================================================================
# sim_T_treatments.R
# Validation of the GENERAL T-treatment MVBCF (order-r truncation) on data with
# known effects. With T=3 and order=2 the engine fits
#   mu, tau_Z1, tau_Z2, tau_Z3, tau_Z1:Z2, tau_Z1:Z3, tau_Z2:Z3.
# The generator has TRUE pairwise interactions {1,2} and {2,3} but NO {1,3};
# a correct model recovers the two real ones and shrinks the spurious {1,3}~0.
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo); library(RcppDist)})
cands <- c("../mvbcfMT", "mvbcfMT", "mvbcf_two_treatments/mvbcfMT")
PKG <- cands[dir.exists(cands)][1]
if (is.na(PKG)) stop("cannot locate mvbcfMT/ ; run from simulations/ or repo root")
dir.create("outputs", showWarnings = FALSE)
sourceCpp(file.path(PKG, "src", "mvbcf_multi_engine.cpp"))
for (f in c("components.R","simulate.R","fit.R","predict.R","bartman.R"))
  source(file.path(PKG, "R", f))

set.seed(7)
sim <- simulate_multi(n = 2000, p = 6, q = 2, Ti = 3, confounded = TRUE, seed = 7)
d   <- sim$data

fit <- fit_mvbcf_multi(d, responses = sim$responses, treatments = sim$treatments,
                       covariates = sim$covariates, order = 2,
                       n_iter = 1200, n_burn = 600, keep_every = 2,
                       n_tree = 80, n_tree_tau = 40, verbose = FALSE)
cat("components:", paste(fit$components$labels, collapse = ", "), "\n\n")

cr <- function(a, b) cor(as.vector(a), as.vector(b))
# main effects
cat("Main effects (corr true vs est):\n")
for (t in 1:3) {
  es <- component_effect(fit, paste0("tau_Z", t), test = FALSE)$mean
  cat(sprintf("  tau_Z%d : %.3f\n", t, cr(sim$truth$main[[t]], es)))
}
# interactions: two real, one spurious
e12 <- component_effect(fit, "tau_Z1:Z2", test = FALSE)$mean
e23 <- component_effect(fit, "tau_Z2:Z3", test = FALSE)$mean
e13 <- component_effect(fit, "tau_Z1:Z3", test = FALSE)$mean
sd_main <- sd(component_effect(fit, "tau_Z2", test = FALSE)$mean)
cat("\nInteractions:\n")
cat(sprintf("  tau_Z1:Z2 (REAL)     corr=%.3f\n", cr(sim$truth$pair[["1:2"]], e12)))
cat(sprintf("  tau_Z2:Z3 (REAL)     corr=%.3f\n", cr(sim$truth$pair[["2:3"]], e23)))
cat(sprintf("  tau_Z1:Z3 (SPURIOUS) sd(est)/sd(main)=%.2f  (near 0 = correctly absent)\n",
            sd(e13) / sd_main))

# joint all-on effect
joint_true <- Reduce(`+`, sim$truth$main) + sim$truth$pair[["1:2"]] + sim$truth$pair[["2:3"]]
joint_est  <- contrast_effect(fit, c(1,1,1), c(0,0,0), test = FALSE)$mean
cat(sprintf("\nJoint (all-on vs all-off) corr=%.3f\n", cr(joint_true, joint_est)))
cat("\nDONE.\n")
