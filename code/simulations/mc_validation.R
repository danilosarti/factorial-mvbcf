# Monte Carlo validation: does the estimator actually capture the effects?
# Over R replications (fresh data each): bias of the ATE, RMSE of the effect
# surface, per-unit CATE 95% coverage, and ATE 95%-interval coverage across reps
# (the clean frequentist calibration check). Plus MCMC ESS on one fit.
suppressMessages({library(Rcpp); library(RcppArmadillo); library(RcppDist)})
PKG <- c("../mvbcfMT","mvbcfMT")[dir.exists(c("../mvbcfMT","mvbcfMT"))][1]
OUT <- "outputs"; dir.create(OUT, showWarnings=FALSE)
sourceCpp(file.path(PKG,"src","mvbcf_multi_engine.cpp"))
for (f in c("components","simulate","fit","predict")) source(file.path(PKG,"R",paste0(f,".R")))

R <- 40; n <- 700; ri <- 1                     # outcome 1
comps <- c("tau_Z1","tau_Z2","tau_Z1:Z2")
truthkey <- list(tau_Z1="m1", tau_Z2="m2", "tau_Z1:Z2"="p12")
csv <- file.path(OUT,"mc_validation_raw.csv")
cat("rep,component,corr,ate_true,ate_est,ate_bias,rmse,cate_cov,ate_covered\n", file=csv)

ess <- function(x){ x<-x-mean(x); n<-length(x); if(sd(x)==0) return(n)
  a<-acf(x, plot=FALSE, lag.max=min(50,n-1))$acf[-1]; k<-which(a<0.05)[1]; if(is.na(k)) k<-length(a)
  n/(1+2*sum(a[1:max(1,k-1)])) }
ess_rec <- NULL

for (r in 1:R) {
  sim <- simulate_multi(n=n, p=5, q=2, Ti=2, seed=100+r)
  fit <- fit_mvbcf_multi(sim$data, responses=sim$responses, treatments=sim$treatments,
                         covariates=sim$covariates, order=2,
                         n_iter=600, n_burn=300, keep_every=2,
                         n_tree=40, n_tree_tau=20, verbose=FALSE)
  truths <- list(m1=sim$truth$main[[1]][,ri], m2=sim$truth$main[[2]][,ri],
                 p12=sim$truth$pair[["1:2"]][,ri])
  for (cl in comps) {
    k  <- match(cl, fit$components$labels)
    cube <- fit$model$predictions[[k]]           # n x q x draws
    tr <- truths[[truthkey[[cl]]]]
    ce <- component_effect(fit, cl, test=FALSE)
    est <- ce$mean[,ri]; lo <- ce$lo[,ri]; hi <- ce$hi[,ri]
    ate_draws <- apply(cube[,ri,], 2, mean)
    aci <- quantile(ate_draws, c(.025,.975)); ate_true <- mean(tr)
    row <- sprintf("%d,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n", r, cl,
      cor(tr,est), ate_true, mean(est), mean(est)-ate_true,
      sqrt(mean((est-tr)^2)), mean(tr>=lo & tr<=hi),
      as.integer(ate_true>=aci[1] & ate_true<=aci[2]))
    cat(row, file=csv, append=TRUE)
    if (r==1) ess_rec <- rbind(ess_rec, data.frame(component=cl, ess_ate=ess(ate_draws), n_draws=length(ate_draws)))
  }
  cat("rep", r, "of", R, "done\n"); flush.console()
}

# aggregate
d <- read.csv(csv)
agg <- do.call(rbind, lapply(split(d, d$component), function(s) data.frame(
  component=s$component[1], reps=nrow(s),
  mean_corr=mean(s$corr),
  mean_ate_bias=mean(s$ate_bias), mc_se_bias=sd(s$ate_bias)/sqrt(nrow(s)),
  rel_bias_pct=100*mean(s$ate_bias)/mean(abs(s$ate_true)),
  mean_rmse=mean(s$rmse),
  cate_coverage=mean(s$cate_cov),
  ate_interval_coverage=mean(s$ate_covered))))
cat("\n===== MONTE CARLO VALIDATION (R=",R,", n=",n,") =====\n", sep="")
print(agg, digits=3, row.names=FALSE)
cat("\nATE-interval coverage should be ~0.95 if calibrated.\n")
cat("\nMCMC ESS of the ATE (one fit):\n"); print(ess_rec, digits=1, row.names=FALSE)
saveRDS(list(raw=d, agg=agg, ess=ess_rec), file.path(OUT,"mc_validation.rds"))
cat("\nALL DONE\n")
