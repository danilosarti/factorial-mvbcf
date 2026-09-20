# =============================================================================
# three_contexts.R  --  the SAME context-agnostic engine on simulated data and
# THREE real domains: clinical trials, agriculture, and economics.
# One modelling call (fit_mvbcf_multi) + one interpretability layer, four datasets.
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo); library(RcppDist); library(ggplot2)
                  library(survival); library(AER); library(MASS)})
PKG <- c("../mvbcfMT","mvbcfMT")[dir.exists(c("../mvbcfMT","mvbcfMT"))][1]
OUT <- "outputs"; dir.create(OUT, showWarnings = FALSE)
sourceCpp(file.path(PKG,"src","mvbcf_multi_engine.cpp"))
for (f in c("components","simulate","fit","predict","bartman","vsup"))
  source(file.path(PKG,"R",paste0(f,".R")))
cr <- function(a,b) cor(as.vector(a), as.vector(b))
summary_rows <- list()
mark <- function(...) { cat(...); cat("\n"); flush.console() }

save_domain <- function(tag, title, fit, resp1, xvar, yvar, benchmark = NULL) {
  ggsave(file.path(OUT, paste0("ctx_",tag,"_forest.pdf")),
         plot_effect_forest(fit, resp1) + labs(subtitle = title), width = 6.5, height = 2.8)
  ggsave(file.path(OUT, paste0("ctx_",tag,"_importance.pdf")),
         plot_component_importance(fit, top = 6), width = 8.5, height = 5.5)
  vv <- tryCatch(plot_effect_vsup_xy(fit, fit$components$labels[2], xvar, yvar, resp1, bins = 7),
                 error = function(e) NULL)
  if (!is.null(vv)) ggsave(file.path(OUT, paste0("ctx_",tag,"_vsup.pdf")), vv, width = 5, height = 4.2)
  # effects table
  tab <- do.call(rbind, lapply(fit$components$labels[!fit$components$is_prognostic], function(l){
    ce <- component_effect(fit, l, test = FALSE); ri <- match(resp1, fit$responses)
    data.frame(effect = l, ate = mean(ce$mean[,ri]), lo = mean(ce$lo[,ri]), hi = mean(ce$hi[,ri]))}))
  cat("\n===", title, "===\n"); print(tab, digits = 3, row.names = FALSE)
  if (!is.null(benchmark)) { cat("benchmark:\n"); print(benchmark) }
  tab$domain <- tag; tab
}

## ---------- 0. SIMULATED (agnostic, known truth) --------------------------
mark("[0] simulated")
set.seed(3)
sim <- simulate_multi(n = 1500, p = 5, q = 2, Ti = 2, seed = 3)
fit0 <- fit_mvbcf_multi(sim$data, responses = sim$responses, treatments = sim$treatments,
                        covariates = sim$covariates, order = 2,
                        n_iter = 900, n_burn = 450, keep_every = 2,
                        n_tree = 60, n_tree_tau = 30, verbose = FALSE)
rc <- c(cr(sim$truth$main[[1]], component_effect(fit0,"tau_Z1",test=FALSE)$mean),
        cr(sim$truth$main[[2]], component_effect(fit0,"tau_Z2",test=FALSE)$mean),
        cr(sim$truth$pair[["1:2"]], component_effect(fit0,"tau_Z1:Z2",test=FALSE)$mean))
cat(sprintf("recovery corr: tau1=%.2f tau2=%.2f tau12=%.2f\n", rc[1],rc[2],rc[3]))
summary_rows$sim <- save_domain("sim","Simulated (known effects)", fit0, "y1", "x1","x2")

## ---------- 1. CLINICAL: rotterdam breast cancer (hormon x chemo) ---------
mark("[1] clinical done-setup")
data(rotterdam, package = "survival")
rt <- rotterdam
set.seed(1); if (nrow(rt) > 1500) rt <- rt[sample(nrow(rt), 1500), ]
rt$logrfs <- log(rt$rtime + 1)          # recurrence-free time (days)
rt$logos  <- log(rt$dtime + 1)          # overall survival time
rt$size_ord <- as.integer(rt$size)       # ordered tumour-size class
cl_cov <- c("age","meno","size_ord","grade","nodes","pgr","er")
fitC <- fit_mvbcf_multi(rt, responses = c("logrfs","logos"),
                        treatments = c("hormon","chemo"), covariates = cl_cov, order = 2,
                        n_iter = 700, n_burn = 350, keep_every = 2,
                        n_tree = 50, n_tree_tau = 25, min_nodesize = 10, verbose = FALSE)
# naive (unadjusted) benchmark for contrast with the propensity-adjusted forest
bmC <- round(coef(lm(logos ~ hormon*chemo, data = rt)), 3)
summary_rows$clin <- save_domain("clinical",
  "Clinical: rotterdam (hormon x chemo), log overall survival", fitC, "logos",
  "age","nodes", benchmark = bmC)

## ---------- 2. AGRICULTURE: npk factorial (N x P x K) ---------------------
mark("[2] agriculture")
data(npk)
np <- npk
np$N <- as.integer(np$N==1); np$P <- as.integer(np$P==1); np$K <- as.integer(np$K==1)
for (b in levels(npk$block)) np[[paste0("blk_",b)]] <- as.integer(npk$block==b)
agcov <- grep("^blk_", names(np), value = TRUE)
fitA <- fit_mvbcf_multi(np, responses = "yield", treatments = c("N","P","K"),
                        covariates = agcov, order = 2,
                        n_iter = 3000, n_burn = 1500, keep_every = 3,
                        n_tree = 40, n_tree_tau = 20, min_nodesize = 2, verbose = FALSE)
bmA <- round(coef(lm(yield ~ N+P+K, data = np)), 2)
summary_rows$ag <- save_domain("agri", "Agriculture: npk (N,P,K), pea yield",
                               fitA, "yield", "blk_2","blk_3", benchmark = bmA)

## ---------- 3. ECONOMICS: hedonic house prices (aircon x prefer) ----------
mark("[3] economics")
data("HousePrices", package = "AER")
hp <- HousePrices
hp$logprice <- log(hp$price); hp$loglot <- log(hp$lotsize)
hp$aircon <- as.integer(hp$aircon == "yes")
hp$prefer <- as.integer(hp$prefer == "yes")
for (v in c("driveway","recreation","fullbase","gasheat"))
  hp[[v]] <- as.integer(hp[[v]] == "yes")
ec_cov <- c("loglot","bedrooms","bathrooms","stories","garage",
            "driveway","recreation","fullbase","gasheat")
fitE <- fit_mvbcf_multi(hp, responses = "logprice",
                        treatments = c("aircon","prefer"), covariates = ec_cov, order = 2,
                        n_iter = 1000, n_burn = 500, keep_every = 2,
                        n_tree = 60, n_tree_tau = 30, verbose = FALSE)
bmE <- round(coef(lm(logprice ~ aircon + prefer, data = hp)), 3)
summary_rows$econ <- save_domain("econ",
  "Economics: hedonic house prices (aircon x prefer), log price",
  fitE, "logprice", "loglot","stories", benchmark = bmE)

## ---------- combined summary ---------------------------------------------
allres <- do.call(rbind, summary_rows)
saveRDS(list(results = allres, sim_recovery = rc),
        file.path(OUT, "three_contexts_results.rds"))
cat("\n==== CROSS-DOMAIN SUMMARY (effect on primary outcome, 95% CrI) ====\n")
print(allres, digits = 3, row.names = FALSE)
cat("\nALL DONE\n")
