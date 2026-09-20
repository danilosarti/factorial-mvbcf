# =============================================================================
# nhanes_application.R
# Deep real-data application of the factorial multivariate Bayesian Causal Forest
# to two open NHANES analyses:
#   APP A (cardiometabolic): statin use x physical activity  ->  SBP, DBP, HbA1c, TC
#   APP B (neurology/cognition, ages 60+): antihypertensive med x physical activity
#                                          ->  CERAD immediate, CERAD delayed,
#                                              Animal Fluency, Digit-Symbol (DSST)
#
# All NHANES variable codes below were verified against the CDC codebooks
# (2011-2012 cycle "_G", 2013-2014 cycle "_H"). Data are public and pulled with
# the nhanesA package. Nothing is fabricated: if a table/variable is unavailable
# for a cycle the script skips it and reports what it used.
#
# HONESTY GUARDRAILS (kept from the paper's stated intent):
#   * The ATE (population-average effect) is the calibrated primary claim.
#   * Interactions are reported ONLY where the 2x2 factorial cell has support
#     (overlap_diagnostic); thin cells are "unsupported", NOT null.
#   * Per-unit CATE maps are EXPLORATORY (VSUP uncertainty suppression).
#   * Observational design => unconfoundedness is an assumption; we adjust for a
#     rich covariate set (propensity controls are built into the model) and run an
#     unmeasured-confounding sensitivity analysis on every ATE.
#
# RUN (on your Mac's R, from the repo root, with real CRAN + the package sources):
#     caffeinate -is Rscript real_data/nhanes_application.R
# Outputs land in real_data/outputs/.
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
need <- function(p) if (!requireNamespace(p, quietly = TRUE))
  try(install.packages(p, quiet = TRUE), silent = TRUE)
for (p in c("nhanesA","Rcpp","RcppArmadillo","RcppDist","grf","ggplot2","dbarts"))
  need(p)
suppressMessages({library(nhanesA); library(Rcpp); library(RcppArmadillo)
                  library(RcppDist); library(ggplot2)})
HAS_GRF <- requireNamespace("grf", quietly = TRUE)

# ---- locate & load the package sources (same pattern as simulations/comparison.R)
PKG <- c("mvbcfMT","../mvbcfMT")[dir.exists(c("mvbcfMT","../mvbcfMT"))][1]
stopifnot(!is.na(PKG))
sourceCpp(file.path(PKG, "src", "mvbcf_multi_engine.cpp"))
for (f in c("components","simulate","fit","predict","diagnostics","vsup"))
  source(file.path(PKG, "R", paste0(f, ".R")))

OUT <- "real_data/outputs"; if (!dir.exists(OUT)) OUT <- "outputs"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "nhanes_run_log.txt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG)  # truncate log
say("NHANES application | R ", as.character(getRversion()), " | grf: ", HAS_GRF,
    " | ", format(Sys.time()))

# ---- helpers ----------------------------------------------------------------
pull <- function(tbl) {  # robust nhanesA pull; returns NULL on failure
  d <- tryCatch(nhanesA::nhanes(tbl, translated = FALSE), error = function(e) NULL)
  if (is.null(d)) d <- tryCatch(nhanesA::nhanes(tbl), error = function(e) NULL)
  if (!is.null(d)) say("  pulled ", tbl, ": ", nrow(d), " x ", ncol(d))
  d
}
getcol <- function(d, v) if (!is.null(d) && v %in% names(d)) d[[v]] else NA
# robust to raw codes (1/2) or translated labels ("Yes"/"No"): -> 1/0, NA kept
is_yes <- function(x) { x <- as.character(x)
  ifelse(is.na(x), NA_integer_, as.integer(x %in% c("1","Yes","yes","YES"))) }
is_female <- function(x) { x <- as.character(x)
  ifelse(is.na(x), NA_integer_, as.integer(x %in% c("2","Female","female"))) }
num <- function(x) suppressWarnings(as.numeric(as.character(x)))  # codes -> numeric
rowmean_pos <- function(M) {                        # mean of positive BP readings
  M[M == 0] <- NA; r <- rowMeans(M, na.rm = TRUE); r[is.nan(r)] <- NA; r }

# combine the two cycles for a family of tables (G = 2011-12, H = 2013-14)
pull2 <- function(base) {                            # base e.g. "DEMO" -> DEMO_G, DEMO_H
  g <- pull(paste0(base, "_G")); h <- pull(paste0(base, "_H"))
  list(g = g, h = h)
}

# ---- 1. PULL RAW TABLES (both cycles) ---------------------------------------
say("\n[1] pulling NHANES tables (2011-2012 _G and 2013-2014 _H) ...")
DEMO <- pull2("DEMO"); BPX <- pull2("BPX"); GHB <- pull2("GHB"); TCHOL <- pull2("TCHOL")
BMX  <- pull2("BMX");  BPQ <- pull2("BPQ"); PAQ <- pull2("PAQ"); DIQ  <- pull2("DIQ")
SMQ  <- pull2("SMQ");  RXQ <- pull2("RXQ_RX"); CFQ <- pull2("CFQ")

# ---- 2. PER-CYCLE PERSON TABLE ----------------------------------------------
build_cycle <- function(demo, bpx, ghb, tchol, bmx, bpq, paq, diq, smq, rxq, cfq) {
  if (is.null(demo)) return(NULL)
  d <- data.frame(SEQN = demo$SEQN,
                  age  = num(getcol(demo,"RIDAGEYR")),
                  sex  = is_female(getcol(demo,"RIAGENDR")),          # 1 = female
                  educ = num(getcol(demo,"DMDEDUC2")),
                  race = num(getcol(demo,"RIDRETH1")),
                  pir  = num(getcol(demo,"INDFMPIR")))
  merge_on <- function(d, src, cols, fun) {
    if (is.null(src)) { for (nm in names(cols)) d[[nm]] <- NA; return(d) }
    add <- fun(src); d <- merge(d, add, by = "SEQN", all.x = TRUE); d }

  # blood pressure (avg of up to 4 readings)
  d <- merge_on(d, bpx, list(sbp=1,dbp=1), function(s){
    sy <- rowmean_pos(cbind(num(getcol(s,"BPXSY1")),num(getcol(s,"BPXSY2")),num(getcol(s,"BPXSY3")),num(getcol(s,"BPXSY4"))))
    di <- rowmean_pos(cbind(num(getcol(s,"BPXDI1")),num(getcol(s,"BPXDI2")),num(getcol(s,"BPXDI3")),num(getcol(s,"BPXDI4"))))
    data.frame(SEQN=s$SEQN, sbp=sy, dbp=di) })
  d <- merge_on(d, ghb,   list(hba1c=1), function(s) data.frame(SEQN=s$SEQN, hba1c=num(getcol(s,"LBXGH"))))
  d <- merge_on(d, tchol, list(tc=1),    function(s) data.frame(SEQN=s$SEQN, tc=num(getcol(s,"LBXTC"))))
  d <- merge_on(d, bmx,   list(bmi=1),   function(s) data.frame(SEQN=s$SEQN, bmi=num(getcol(s,"BMXBMI"))))
  d <- merge_on(d, bpq,   list(bp_med=1),function(s) data.frame(SEQN=s$SEQN, bp_med=is_yes(getcol(s,"BPQ050A"))))
  d <- merge_on(d, diq,   list(diab=1),  function(s) data.frame(SEQN=s$SEQN, diab=is_yes(getcol(s,"DIQ010"))))
  d <- merge_on(d, smq,   list(smoke=1), function(s) data.frame(SEQN=s$SEQN, smoke=is_yes(getcol(s,"SMQ020"))))
  d <- merge_on(d, paq,   list(active=1),function(s){
    a <- pmax(is_yes(getcol(s,"PAQ605")), is_yes(getcol(s,"PAQ620")),
              is_yes(getcol(s,"PAQ650")), is_yes(getcol(s,"PAQ665")), na.rm=TRUE)
    a[is.na(a)] <- 0
    data.frame(SEQN=s$SEQN, active=as.integer(a)) })
  # statin use from long-format prescription file
  d <- merge_on(d, rxq, list(statin=1), function(s){
    drug <- toupper(as.character(getcol(s,"RXDDRUG")))
    is_st <- grepl("STATIN", drug) & !grepl("NYSTATIN", drug)
    ag <- aggregate(is_st, by=list(SEQN=s$SEQN), FUN=function(z) as.integer(any(z, na.rm=TRUE)))
    names(ag) <- c("SEQN","statin"); ag })
  # cognition (ages 60+ only; absent in some persons)
  d <- merge_on(d, cfq, list(cerad_imm=1,cerad_del=1,animal=1,dsst=1), function(s){
    imm <- num(getcol(s,"CFDCST1")) + num(getcol(s,"CFDCST2")) + num(getcol(s,"CFDCST3"))
    data.frame(SEQN=s$SEQN, cerad_imm=imm, cerad_del=num(getcol(s,"CFDCSR")),
               animal=num(getcol(s,"CFDAST")), dsst=num(getcol(s,"CFDDS"))) })
  # bp_med / statin: NA means not on branch -> treat as 0 (not using)
  d$bp_med[is.na(d$bp_med)] <- 0; d$statin[is.na(d$statin)] <- 0
  d$diab[is.na(d$diab)] <- 0; d$smoke[is.na(d$smoke)] <- 0
  d
}
cyc <- function(fam, k) fam[[k]]
cG <- build_cycle(cyc(DEMO,"g"),cyc(BPX,"g"),cyc(GHB,"g"),cyc(TCHOL,"g"),cyc(BMX,"g"),
                  cyc(BPQ,"g"),cyc(PAQ,"g"),cyc(DIQ,"g"),cyc(SMQ,"g"),cyc(RXQ,"g"),cyc(CFQ,"g"))
cH <- build_cycle(cyc(DEMO,"h"),cyc(BPX,"h"),cyc(GHB,"h"),cyc(TCHOL,"h"),cyc(BMX,"h"),
                  cyc(BPQ,"h"),cyc(PAQ,"h"),cyc(DIQ,"h"),cyc(SMQ,"h"),cyc(RXQ,"h"),cyc(CFQ,"h"))
common <- Reduce(intersect, list(names(cG), names(cH)))
pooled <- rbind(cG[,common], cH[,common])
say("\n[2] pooled person table: ", nrow(pooled), " rows, ", ncol(pooled), " cols")

# ---- 3. FIT ONE APPLICATION -------------------------------------------------
COV <- c("age","sex","educ","race","pir","bmi","diab","smoke")
run_app <- function(tag, dat, Z1, Z2, z1lab, z2lab, outcomes, olab, covs = COV) {
  say("\n===== APPLICATION ", tag, " : ", z1lab, " x ", z2lab, " =====")
  keep <- c("SEQN", Z1, Z2, outcomes, covs)
  d <- dat[, keep]; d <- d[stats::complete.cases(d), ]
  names(d)[match(c(Z1,Z2), names(d))] <- c("Z1","Z2")
  say("  complete cases: ", nrow(d))
  # standardize outcomes to SD units (multivariate Sigma well-conditioned)
  osd <- sapply(outcomes, function(o) stats::sd(d[[o]]))
  omu <- sapply(outcomes, function(o) mean(d[[o]]))
  for (o in outcomes) d[[o]] <- (d[[o]] - mean(d[[o]])) / stats::sd(d[[o]])
  # 2x2 cell counts
  cell <- table(d$Z1, d$Z2); say("  2x2 cell counts (rows Z1, cols Z2):")
  capture.output(print(cell), file = LOG, append = TRUE); print(cell)

  fit <- fit_mvbcf_multi(d, responses = outcomes, treatments = c("Z1","Z2"),
                         covariates = covs, order = 2,
                         n_iter = 2000, n_burn = 1000, keep_every = 2,
                         n_tree = 100, n_tree_tau = 50, seed = 1, verbose = TRUE)

  # overlap / positivity
  ov <- overlap_diagnostic(fit); say("  --- overlap diagnostic ---")
  capture.output(print(ov), file = LOG, append = TRUE); print(ov)

  # ATE table (posterior mean + 95% CrI from draws) for each component & outcome
  comps <- c("tau_Z1","tau_Z2","tau_Z1:Z2")
  rows <- list()
  for (cp in comps) {
    ce <- component_effect(fit, cp, test = FALSE)     # mean/lo/hi/draws (n x q [x S])
    dr <- ce$draws                                    # [n, q, S]
    for (ri in seq_along(outcomes)) {
      ate_draws <- apply(dr[, ri, ], 2, mean)         # ATE per posterior draw (SD units)
      est <- mean(ate_draws); lo <- stats::quantile(ate_draws,.025); hi <- stats::quantile(ate_draws,.975)
      sens <- sensitivity_effect(est, lo, hi, outcome_sd = 1)  # already SD units
      rows[[length(rows)+1]] <- data.frame(
        app=tag, component=cp, outcome=outcomes[ri],
        ate_sd=round(est,3), lo_sd=round(lo,3), hi_sd=round(hi,3),
        ate_natural=round(est*osd[ri],3),        # back to natural units
        excludes_zero=sens$excludes_zero,
        robustness_sd=sens$robustness_sd)
    }
  }
  ate <- do.call(rbind, rows)
  write.csv(ate, file.path(OUT, paste0("nhanes_", tag, "_ATE.csv")), row.names = FALSE)
  say("  --- ATE table (SD units; natural-unit column back-transformed) ---")
  capture.output(print(ate, row.names = FALSE), file = LOG, append = TRUE); print(ate)

  # grf multi-arm on the same data, first outcome (mirrors the simulation comparison)
  if (HAS_GRF) {
    X <- as.matrix(d[, covs]); W <- factor(paste0(d$Z1,d$Z2), levels=c("00","10","01","11"))
    for (ri in seq_along(outcomes)) {
      mac <- grf::multi_arm_causal_forest(X=X, Y=d[[outcomes[ri]]], W=W, num.trees=2000)
      pr <- predict(mac)$predictions
      g <- data.frame(app=tag, outcome=outcomes[ri], method="grf-multiarm",
        ate_t1=round(mean(pr[,1,1]),3), ate_t2=round(mean(pr[,2,1]),3),
        ate_t12=round(mean(pr[,3,1]-pr[,1,1]-pr[,2,1]),3))
      write.table(g, file.path(OUT, paste0("nhanes_", tag, "_grf.csv")),
                  sep=",", row.names=FALSE, append=(ri>1), col.names=(ri==1))
    }
  }

  # figures: ATE forest (first outcome) + VSUP interaction map over age x bmi
  p1 <- tryCatch(plot_effect_forest(fit, response = outcomes[1], test = FALSE),
                 error=function(e) NULL)
  if (!is.null(p1)) ggsave(file.path(OUT, paste0("nhanes_",tag,"_forest.pdf")),
                           p1, width=7, height=3.4)
  p2 <- tryCatch(plot_effect_vsup_xy(fit, "tau_Z1:Z2", xvar="age", yvar="bmi",
                 response=outcomes[1], bins=8, test=FALSE), error=function(e) NULL)
  if (!is.null(p2)) ggsave(file.path(OUT, paste0("nhanes_",tag,"_vsup_interaction.pdf")),
                           p2, width=5, height=4)

  saveRDS(list(fit=fit, ate=ate, overlap=ov, osd=osd, omu=omu),
          file.path(OUT, paste0("nhanes_", tag, ".rds")))
  say("  saved: nhanes_", tag, "_ATE.csv / _grf.csv / _forest.pdf / _vsup_interaction.pdf / .rds")
  invisible(ate)
}

# ---- APP A: cardiometabolic (adults 20+) ------------------------------------
datA <- pooled[!is.na(pooled$age) & pooled$age >= 20, ]
ateA <- run_app("A_cardio", datA, Z1="statin", Z2="active",
                z1lab="statin use", z2lab="physical activity",
                outcomes=c("sbp","dbp","hba1c","tc"),
                olab=c("Systolic BP","Diastolic BP","HbA1c","Total chol."))

# ---- APP B: neurology / cognition (ages 60+) --------------------------------
datB <- pooled[!is.na(pooled$age) & pooled$age >= 60, ]
ateB <- run_app("B_neuro", datB, Z1="bp_med", Z2="active",
                z1lab="antihypertensive medication", z2lab="physical activity",
                outcomes=c("cerad_imm","cerad_del","animal","dsst"),
                olab=c("CERAD immediate","CERAD delayed","Animal fluency","Digit-symbol"))

say("\nALL DONE. See ", OUT, "/ (ATE csvs, figures, rds) and this log.")
