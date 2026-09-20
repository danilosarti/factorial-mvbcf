# =============================================================================
# accord_application.R  --  ITEM 3: randomized 2x2 factorial confirmation
# Factorial MVBCF on the ACCORD trial (NHLBI BioLINCC), the strongest possible
# identification of a treatment INTERACTION because the crossing is RANDOMIZED.
#
# STATUS: TEMPLATE, pending data. ACCORD is CONTROLLED ACCESS: it requires an
# approved BioLINCC Data Use Agreement (weeks of lead time). This script cannot
# be run until the data are in hand; the variable names marked  <<< SET >>>  must
# be matched to the delivered data dictionary (BioLINCC ships SAS + a codebook).
# WHY randomized matters: with ACCORD the interaction is identified WITHOUT the
# unconfoundedness assumption the NHANES analysis leans on -- it removes the #1
# caveat a referee raises against the observational application.
#
# DESIGN (ACCORD): every participant randomized to intensive vs standard GLYCEMIA,
# and, in the BP cohort, additionally to intensive vs standard BLOOD PRESSURE ->
# a clean randomized 2x2 (glycemia x BP). (The lipid cohort gives glycemia x
# fenofibrate as an alternative 2x2.) We estimate the two main effects and their
# interaction on a multivariate panel of continuous outcomes, jointly.
#
# RUN (once the DUA is approved and the CSVs are local, on your Mac's R):
#     caffeinate -is Rscript real_data/accord_application.R
# =============================================================================

options(stringsAsFactors = FALSE)
suppressMessages({library(Rcpp); library(RcppArmadillo); library(RcppDist)})
HAS_GRF <- requireNamespace("grf", quietly = TRUE)

PKG <- c("mvbcfMT","../mvbcfMT")[dir.exists(c("mvbcfMT","../mvbcfMT"))][1]
stopifnot(!is.na(PKG))
sourceCpp(file.path(PKG,"src","mvbcf_multi_engine.cpp"))
for (f in c("components","simulate","fit","predict","diagnostics","vsup"))
  source(file.path(PKG,"R",paste0(f,".R")))
OUT <- "real_data/outputs"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# ---- 1. LOAD ACCORD (edit paths to the BioLINCC delivery) -------------------
#   BioLINCC ships several tables; you typically need the randomization/treatment
#   assignment file, the baseline covariate file, and a follow-up measurements
#   file. Merge them on the participant id (MaskID in ACCORD).
ACCORD_DIR <- "~/data/accord"                                   # <<< SET >>>
rand  <- read.csv(file.path(ACCORD_DIR, "accord_key.csv"))      # <<< SET >>> assignment
base  <- read.csv(file.path(ACCORD_DIR, "bloodpressure.csv"))   # <<< SET >>> baseline+labs
labs  <- read.csv(file.path(ACCORD_DIR, "lipids.csv"))          # <<< SET >>> labs
# ... merge as needed into one person-level frame `d0` keyed by MaskID.

# ---- 2. BUILD THE 2x2 AND THE OUTCOME PANEL ---------------------------------
# Map to the analysis frame. Treatments are RANDOMIZED indicators (0/1):
#   Z1 = intensive glycemia   (from the glycemia arm variable)   <<< SET >>>
#   Z2 = intensive BP         (from the BP arm variable)         <<< SET >>>
# Outcomes: a panel of correlated continuous endpoints at a fixed follow-up
# (e.g., 12 months), standardized before fitting:
#   hba1c, sbp, ldl, bmi   (choose 3-4 available, correlated)    <<< SET >>>
# Covariates: baseline versions of the outcomes + age, sex, race, smoking,
# CVD history, baseline eGFR, etc.                                <<< SET >>>
#
# d <- data.frame(MaskID=..., Z1=..., Z2=...,
#                 hba1c=..., sbp=..., ldl=..., bmi=...,
#                 age=..., sex=..., ... )
# d <- d[complete.cases(d), ]
# outcomes <- c("hba1c","sbp","ldl","bmi")
# covs     <- c("age","sex","race","smoke","cvd_hx","egfr",
#               paste0("base_", outcomes))
# osd <- sapply(outcomes, function(o) sd(d[[o]]))
# for (o in outcomes) d[[o]] <- scale(d[[o]])[,1]

# ---- 3. FIT (identical call to the NHANES pipeline) -------------------------
# fit <- fit_mvbcf_multi(d, responses=outcomes, treatments=c("Z1","Z2"),
#                        covariates=covs, order=2,
#                        n_iter=2000, n_burn=1000, keep_every=2,
#                        n_tree=100, n_tree_tau=50, seed=1, verbose=TRUE)
#
# Because assignment is randomized, the propensity controls are ~0.5 by design;
# the overlap diagnostic should show all four cells well populated and NO
# near-violations -- the point of using a factorial RCT.
# ov <- overlap_diagnostic(fit); print(ov)
#
# ATE table with 95% CrI from the posterior draws (same code as NHANES):
# for (cp in c("tau_Z1","tau_Z2","tau_Z1:Z2")) {
#   ce <- component_effect(fit, cp, test=FALSE); dr <- ce$draws
#   for (ri in seq_along(outcomes)) {
#     ate <- apply(dr[,ri,],2,mean)
#     cat(cp, outcomes[ri], round(mean(ate),3),
#         round(quantile(ate,.025),3), round(quantile(ate,.975),3), "\n") } }
#
# Head-to-head vs grf multi-arm on the same randomized data, and the VSUP
# interaction map, exactly as in nhanes_application.R.

# ---- 4. THE PAYOFF ----------------------------------------------------------
# Randomization identifies tau_Z1:Z2 without unconfoundedness: a supported,
# nonzero glycemia x BP interaction on the multivariate endpoint panel would be
# a genuinely causal factorial-interaction finding -- the strongest possible
# version of the paper's application, and the one that most raises the ceiling
# for AOAS. A null interaction is equally publishable (a calibrated, randomized
# "no synergy" statement on correlated endpoints).

cat("TEMPLATE ONLY -- fill the  <<< SET >>>  mappings after the BioLINCC DUA,",
    "then this runs identically to nhanes_application.R.\n")
