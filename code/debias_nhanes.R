# =============================================================================
# debias_nhanes.R -- apply the cross-fitted debiased estimator to the saved
# NHANES fits (no refitting: loads the .rds produced by nhanes_application*.R).
# Produces the debiased average effects + efficient SE + Wald CI next to the raw
# posterior plug-in, for both applications. Run on your Mac's R.
#   caffeinate -is Rscript real_data/debias_nhanes.R
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo); library(RcppDist)})
PKG <- c("mvbcfMT","../mvbcfMT")[dir.exists(c("mvbcfMT","../mvbcfMT"))][1]
sourceCpp(file.path(PKG,"src","mvbcf_multi_engine.cpp"))
for (f in c("components","simulate","fit","predict","diagnostics","vsup","debias"))
  source(file.path(PKG,"R",paste0(f,".R")))     # <- put debias.R in mvbcfMT/R/
OUT <- "real_data/outputs"; if (!dir.exists(OUT)) OUT <- "outputs"

for (tag in c("A_cardio","B_neuro")) {
  rds <- file.path(OUT, paste0("nhanes_", tag, ".rds"))
  if (!file.exists(rds)) { cat("missing", rds, "\n"); next }
  fit <- readRDS(rds)$fit
  cat("\n===== debiased ATEs:", tag, "=====\n")
  db <- debiased_ate(fit, folds = 5, propensity = "auto", level = 0.95, seed = 1)
  print(db)
  write.csv(db, file.path(OUT, paste0("nhanes_", tag, "_debiased.csv")), row.names = FALSE)
}
cat("\nDONE. Compare *_debiased.csv (debiased + efficient SE) with *_ATE.csv (raw posterior).\n")
