# =============================================================================
# sim_two_treatments.R
# Validation of the TWO-treatment factorial MVBCF on data with KNOWN effects.
# Fits (mu, tau1, tau2, tau12) and checks recovery of each surface, the ATEs,
# the interaction, and the credible-interval coverage.
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo); library(RcppDist)})
have_gg <- requireNamespace("ggplot2", quietly = TRUE)

# Resolve the package folder whether this is run from simulations/ or repo root.
cands <- c("../mvbcfMT", "mvbcfMT", "mvbcf_two_treatments/mvbcfMT")
PKG <- cands[dir.exists(cands)][1]
if (is.na(PKG)) stop("cannot locate mvbcfMT/ ; run from simulations/ or repo root")
dir.create("outputs", showWarnings = FALSE)
sourceCpp(file.path(PKG, "src", "mvbcf_multi_engine.cpp"))
for (f in c("components.R","simulate.R","fit.R","predict.R","bartman.R"))
  source(file.path(PKG, "R", f))

set.seed(11)
sim <- simulate_met2(n_gen = 40, n_env = 10, confounded = TRUE, seed = 11)
d   <- sim$data
cat(sprintf("rows=%d  Z1=%.2f  Z2=%.2f  both-on=%.2f\n",
            nrow(d), mean(d$Z1), mean(d$Z2), mean(d$Z1==1 & d$Z2==1)))

fit <- fit_mvbcf2(d, responses = c("yield_kg_ha","tgw_g","plant_height_cm"),
                  treatments = c("Z1","Z2"), gen = "gen", env = "env",
                  n_iter = 1500, n_burn = 750, keep_every = 2,
                  n_tree = 100, n_tree_tau = 50, seed = 1, verbose = FALSE)
print(fit)

# --- recovery on training rows (ground truth available) ----------------------
tru <- function(pfx) as.matrix(d[, paste0("true_", pfx, c("_yield","_tgw","_height"))])
est <- function(lab) component_effect(fit, lab, test = FALSE)$mean
tab <- data.frame(
  component = c("tau1","tau2","tau12"),
  corr = c(cor(c(tru("tau1")), c(est("tau_Z1"))),
           cor(c(tru("tau2")), c(est("tau_Z2"))),
           cor(c(tru("tau12")),c(est("tau_Z1:Z2")))),
  ate_true = c(mean(tru("tau1")[,1]), mean(tru("tau2")[,1]), mean(tru("tau12")[,1])),
  ate_est  = c(mean(est("tau_Z1")[,1]), mean(est("tau_Z2")[,1]), mean(est("tau_Z1:Z2")[,1])))
cat("\n--- yield-scale recovery ---\n"); print(tab, digits = 3, row.names = FALSE)

# --- 95% credible-interval coverage for tau1 (yield) -------------------------
c1 <- component_effect(fit, "tau_Z1", test = FALSE)
cov <- mean(tru("tau1")[,1] >= c1$lo[,1] & tru("tau1")[,1] <= c1$hi[,1])
cat(sprintf("\n95%% CI coverage, tau1(yield): %.2f\n", cov))

# --- counterfactual G x E grid incl. never-grown cells -----------------------
pe <- predict_effects(fit, test = TRUE)
cat(sprintf("grid cells=%d  never-grown=%d\n", nrow(pe), sum(pe$observed==0)))
write.csv(pe, "outputs/two_treatments_gxe_effects.csv", row.names = FALSE)

# --- optional plots ----------------------------------------------------------
if (have_gg) {
  library(ggplot2)
  pd <- rbind(
    data.frame(comp="tau1 (drought)",  true=tru("tau1")[,1],  est=est("tau_Z1")[,1]),
    data.frame(comp="tau2 (inoculant)",true=tru("tau2")[,1],  est=est("tau_Z2")[,1]),
    data.frame(comp="tau12 (synergy)", true=tru("tau12")[,1], est=est("tau_Z1:Z2")[,1]))
  g <- ggplot(pd, aes(true, est)) + geom_point(alpha=.4, colour="#2166AC") +
    geom_abline(slope=1, intercept=0, colour="grey40", linetype=2) +
    facet_wrap(~comp, scales="free") +
    labs(title="Two-treatment MVBCF: true vs estimated effect (yield)",
         x="true effect", y="posterior-mean effect") + theme_minimal(base_size=10)
  ggsave("outputs/two_treatments_recovery.pdf", g, width=9, height=3.2)
  cat("wrote outputs/two_treatments_recovery.pdf\n")
}
cat("\nDONE.\n")
