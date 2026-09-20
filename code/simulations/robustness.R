# =============================================================================
# robustness.R  --  closing the asterisk:
#  (A) sample-size sweep: does bias -> 0 and coverage -> nominal as n grows?
#  (B) misspecification stress tests: strong confounding, an UNMEASURED
#      confounder, heavy-tailed + heteroscedastic errors, and under-specified
#      interaction order (fit order 1 while the truth has a pairwise interaction).
# Parallelised over cores; each replication writes its own part-file (crash-safe,
# monitorable). Metrics vs known truth: ATE relative bias, RMSE, per-unit (CATE)
# 95% coverage, and ATE-interval 95% coverage.
# =============================================================================
suppressMessages({library(Rcpp); library(RcppArmadillo); library(RcppDist); library(parallel)})
PKG <- c("../mvbcfMT","mvbcfMT")[dir.exists(c("../mvbcfMT","mvbcfMT"))][1]
OUT <- "outputs"; PARTS <- file.path(OUT,"rob_parts"); dir.create(PARTS, recursive=TRUE, showWarnings=FALSE)
sourceCpp(file.path(PKG,"src","mvbcf_multi_engine.cpp"))
for (f in c("components","simulate","fit","predict")) source(file.path(PKG,"R",paste0(f,".R")))

## ---- data generator with misspecification levers -------------------------
gen <- function(scenario, n, seed) {
  set.seed(seed)
  p <- 5; X <- matrix(rnorm(n*p), n, p); colnames(X) <- paste0("x",1:p)
  x1<-X[,1]; x2<-X[,2]; x3<-X[,3]; u <- rnorm(n)      # u = potential unmeasured confounder
  Ti <- if (grepl("^T3", scenario)) 3 else 2
  strong <- scenario=="strong_confounding"; uc <- if (scenario=="unmeasured") 1.8 else 0
  a1 <- if(strong) 2.2 else 0.6; a3 <- if(strong) 1.6 else 0.4
  e1 <- plogis(a1*x1 + a3*x3 + uc*u - 0.1);  Z1 <- rbinom(n,1,e1)
  e2 <- plogis(-0.9*x2 + 0.5*x1 + 0.8*uc*u - 0.1); Z2 <- rbinom(n,1,e2)
  Z3 <- NULL; if (Ti==3){ e3 <- plogis(0.5*x3 - 0.2); Z3 <- rbinom(n,1,e3) }
  mu <- 2 + 1.5*x1 - 0.4*x2^2 + 0.8*x3 + (if(scenario=="unmeasured") 2.0*u else 0)
  t1 <- 1.5 + 1.2*x1*(x1>0) + 0.9*(x2>0); t2 <- 1.0 - 0.8*x2 + 0.6*x3
  t12 <- 0.9*(x1>0)*(x2>0); t3 <- 1.2 + 0.7*x3; t23 <- -0.8*x2
  q <- 2; qs <- c(1, 1.6)
  Ey <- sapply(1:q, function(k){ s<-qs[k]
    v <- s*mu + Z1*s*t1 + Z2*s*t2 + Z1*Z2*s*t12
    if (Ti==3) v <- v + Z3*s*t3 + Z2*Z3*s*t23     # real 2:3 pairwise; NO 3-way
    v })
  Sig <- matrix(c(1,0.4,0.4,1),2)
  if (scenario=="heavytail") {
    E <- (matrix(rt(n*q, df=3)/sqrt(3), n) %*% chol(Sig)) * (1 + 0.8*abs(x1))
  } else E <- matrix(rnorm(n*q), n) %*% chol(Sig)
  y <- Ey + E; colnames(y) <- paste0("y",1:q)
  df <- data.frame(y, Z1=Z1, Z2=Z2); if (Ti==3) df$Z3 <- Z3; df <- cbind(df, X)
  truth <- list("tau_Z1"=outer(t1,qs), "tau_Z2"=outer(t2,qs), "tau_Z1:Z2"=outer(t12,qs))
  if (Ti==3){ truth[["tau_Z3"]]<-outer(t3,qs); truth[["tau_Z2:Z3"]]<-outer(t23,qs) }
  list(data=df, responses=colnames(y),
       treatments=if(Ti==3) c("Z1","Z2","Z3") else c("Z1","Z2"),
       covariates=colnames(X), truth=truth,
       overlap=min(pmin(e1,1-e1)))          # worst-case propensity overlap
}

metrics <- function(fit, truth, ri=1) {
  do.call(rbind, lapply(intersect(names(truth), fit$components$labels), function(cl){
    k <- match(cl, fit$components$labels); cube <- fit$model$predictions[[k]]
    tr <- truth[[cl]][,ri]; ce <- component_effect(fit, cl, test=FALSE)
    est<-ce$mean[,ri]; lo<-ce$lo[,ri]; hi<-ce$hi[,ri]
    ad <- apply(cube[,ri,],2,mean); aci <- quantile(ad,c(.025,.975)); at <- mean(tr)
    data.frame(component=cl, corr=cor(tr,est), ate_true=at, ate_est=mean(est),
      ate_bias=mean(est)-at, rmse=sqrt(mean((est-tr)^2)),
      cate_cov=mean(tr>=lo & tr<=hi),
      ate_covered=as.integer(at>=aci[1] & at<=aci[2]))
  }))
}

run_one <- function(spec) {
  pf <- file.path(PARTS, sprintf("%s_n%d_r%02d.csv", spec$scenario, spec$n, spec$rep))
  if (file.exists(pf)) return(invisible(NULL))
  out <- tryCatch({
    g <- gen(spec$scenario, spec$n, seed=1000*spec$rep + spec$n)
    fit <- fit_mvbcf_multi(g$data, responses=g$responses, treatments=g$treatments,
             covariates=g$covariates, order=spec$order,
             n_iter=500, n_burn=250, keep_every=2, n_tree=35, n_tree_tau=18,
             min_nodesize=5, verbose=FALSE)
    m <- metrics(fit, g$truth)
    m$scenario<-spec$scenario; m$n<-spec$n; m$rep<-spec$rep; m$order<-spec$order; m$overlap<-g$overlap
    m
  }, error=function(e) data.frame(component=NA, scenario=spec$scenario, n=spec$n, rep=spec$rep, err=conditionMessage(e)))
  write.csv(out, pf, row.names=FALSE)
  invisible(NULL)
}

## ---- build the job list --------------------------------------------------
specs <- list()
# (A) sample-size sweep (well-specified, T=2, order 2)
for (nn in c(300,700,1500,3000)) for (r in 1:12)
  specs[[length(specs)+1]] <- list(scenario="wellspec", n=nn, rep=r, order=2)
# (B) misspecification at n=1200, order 2
for (sc in c("strong_confounding","unmeasured","heavytail")) for (r in 1:12)
  specs[[length(specs)+1]] <- list(scenario=sc, n=1200, rep=r, order=2)
# (C) interaction-order underspecification on a T=3 truth (real 2:3 pairwise)
for (r in 1:12) specs[[length(specs)+1]] <- list(scenario="T3_correct", n=1200, rep=r, order=2)
for (r in 1:12) specs[[length(specs)+1]] <- list(scenario="T3_underspec", n=1200, rep=r, order=1)

cat("total jobs:", length(specs), " on", detectCores(), "cores\n"); flush.console()
invisible(mclapply(specs, run_one, mc.cores=max(1, detectCores()), mc.preschedule=FALSE))

## ---- aggregate -----------------------------------------------------------
parts <- list.files(PARTS, pattern="\\.csv$", full.names=TRUE)
d <- do.call(rbind, lapply(parts, function(f) tryCatch(read.csv(f), error=function(e) NULL)))
d <- d[!is.na(d$component), ]
key <- interaction(d$scenario, d$n, d$component, drop=TRUE)
agg <- do.call(rbind, lapply(split(d, key), function(s) data.frame(
  scenario=s$scenario[1], n=s$n[1], component=s$component[1], reps=nrow(s),
  corr=round(mean(s$corr),3),
  rel_bias_pct=round(100*mean(s$ate_bias)/mean(abs(s$ate_true)+1e-9),1),
  rmse=round(mean(s$rmse),3),
  cate_cov=round(mean(s$cate_cov),3),
  ate_cov=round(mean(s$ate_covered),3),
  overlap=round(mean(s$overlap),3))))
agg <- agg[order(agg$scenario, agg$n, agg$component), ]
saveRDS(list(raw=d, agg=agg), file.path(OUT,"robustness.rds"))
write.csv(agg, file.path(OUT,"robustness_summary.csv"), row.names=FALSE)
cat("\n===== ROBUSTNESS SUMMARY =====\n"); print(agg, row.names=FALSE)
cat("\nALL DONE\n")
