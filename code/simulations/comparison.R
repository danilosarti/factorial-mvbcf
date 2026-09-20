# =============================================================================
# comparison.R -- head-to-head vs the natural alternatives, on a common DGP with
# KNOWN heterogeneous main + interaction effects and correlated outcomes.
# Run on a machine with grf/glmnet/randomForest installed (real CRAN access).
# Competitors:
#   (1) MVBCF-factorial (ours): mu, tau1, tau2, tau12, joint multivariate.
#   (1b) ablation: ours per-outcome (diagonal Sigma).
#   (2) Separate single-treatment causal forests (ours, order 1, one per Z).
#   (3) grf::multi_arm_causal_forest  -- REAL generalized random forest.
#   (4) Random-forest 4-arm T-learner (naive nonparametric baseline).
#   (5) Regularized linear factorial (glmnet); FactorHet best-effort if present.
# Metric vs truth: RMSE and correlation of the per-unit effect; ATE bias.
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
need <- function(p) if(!requireNamespace(p, quietly=TRUE))
  try(install.packages(p, quiet=TRUE), silent=TRUE)
for (p in c("Rcpp","RcppArmadillo","RcppDist","glmnet","randomForest","grf","ggplot2")) need(p)

suppressMessages({library(Rcpp);library(RcppArmadillo);library(RcppDist)
                  library(parallel);library(randomForest);library(glmnet)})
HAS_GRF <- requireNamespace("grf", quietly=TRUE)
HAS_FH  <- requireNamespace("FactorHet", quietly=TRUE)
cat("R", as.character(getRversion()),
    "| grf:", HAS_GRF, "| FactorHet:", HAS_FH, "| glmnet: TRUE\n")

PKG<-c("../mvbcfMT","mvbcfMT")[dir.exists(c("../mvbcfMT","mvbcfMT"))][1]
OUT<-"outputs"; PARTS<-file.path(OUT,"cmp_parts"); dir.create(PARTS,recursive=TRUE,showWarnings=FALSE)
sourceCpp(file.path(PKG,"src","mvbcf_multi_engine.cpp"))
for(f in c("components","simulate","fit","predict")) source(file.path(PKG,"R",paste0(f,".R")))

rmse<-function(a,b)sqrt(mean((a-b)^2)); crr<-function(a,b)suppressWarnings(cor(a,b))
anova4<-function(m00,m10,m01,m11) list(t1=m10-m00,t2=m01-m00,t12=m11-m10-m01+m00)
gdesign<-function(X,z1,z2){ z12<-z1*z2; cbind(X, Z1=z1,Z2=z2,Z12=z12, X*z1, X*z2, X*z12) }

run_rep<-function(r){
  pf<-file.path(PARTS,sprintf("rep%02d.csv",r)); if(file.exists(pf))return(invisible())
  out<-tryCatch({
   sim<-simulate_multi(n=1200,p=5,q=2,Ti=2,seed=500+r); d<-sim$data
   X<-as.matrix(d[,sim$covariates]); Z1<-d$Z1; Z2<-d$Z2; ri<-1
   y<-d[[sim$responses[ri]]]
   tr<-list(t1=sim$truth$main[[1]][,ri],t2=sim$truth$main[[2]][,ri],t12=sim$truth$pair[["1:2"]][,ri])
   rows<-list(); add<-function(est,comp,e){rows[[length(rows)+1]]<<-data.frame(
      rep=r,estimator=est,component=comp,rmse=rmse(e,tr[[comp]]),corr=crr(e,tr[[comp]]),
      ate_bias=mean(e)-mean(tr[[comp]]))}

   ## (1) ours factorial (joint q=2)
   f<-fit_mvbcf_multi(d,responses=sim$responses,treatments=c("Z1","Z2"),covariates=sim$covariates,
        order=2,n_iter=500,n_burn=250,keep_every=2,n_tree=35,n_tree_tau=18,verbose=FALSE)
   add("MVBCF-factorial","t1",component_effect(f,"tau_Z1",test=FALSE)$mean[,ri])
   add("MVBCF-factorial","t2",component_effect(f,"tau_Z2",test=FALSE)$mean[,ri])
   add("MVBCF-factorial","t12",component_effect(f,"tau_Z1:Z2",test=FALSE)$mean[,ri])

   ## (1b) ablation: ours per-outcome (q=1)
   f1<-fit_mvbcf_multi(d,responses=sim$responses[ri],treatments=c("Z1","Z2"),
        covariates=sim$covariates,order=2,n_iter=500,n_burn=250,keep_every=2,
        n_tree=35,n_tree_tau=18,verbose=FALSE)
   add("MVBCF-per-outcome","t1",component_effect(f1,"tau_Z1",test=FALSE)$mean[,1])
   add("MVBCF-per-outcome","t12",component_effect(f1,"tau_Z1:Z2",test=FALSE)$mean[,1])

   ## (2) separate single-treatment forests (order 1 each)
   fa<-fit_mvbcf_multi(d,responses=sim$responses,treatments=c("Z1"),covariates=sim$covariates,
        order=1,n_iter=500,n_burn=250,keep_every=2,n_tree=35,n_tree_tau=18,verbose=FALSE)
   fb<-fit_mvbcf_multi(d,responses=sim$responses,treatments=c("Z2"),covariates=sim$covariates,
        order=1,n_iter=500,n_burn=250,keep_every=2,n_tree=35,n_tree_tau=18,verbose=FALSE)
   add("Separate-forests","t1",component_effect(fa,"tau_Z1",test=FALSE)$mean[,ri])
   add("Separate-forests","t2",component_effect(fb,"tau_Z2",test=FALSE)$mean[,ri])
   add("Separate-forests","t12",rep(0,nrow(d)))  # interaction not estimable -> zero predictor

   ## (3) grf::multi_arm_causal_forest -- REAL grf. Arms = 4 cells of the 2x2.
   if (HAS_GRF) {
     W<-factor(paste0(Z1,Z2), levels=c("00","10","01","11"))
     mac<-grf::multi_arm_causal_forest(X=X, Y=y, W=W, num.trees=2000)
     pr<-predict(mac)$predictions   # [n, contrasts(=10,01,11 vs 00), outcomes(=1)]
     c10<-pr[,1,1]; c01<-pr[,2,1]; c11<-pr[,3,1]
     add("grf-multiarm","t1",c10); add("grf-multiarm","t2",c01)
     add("grf-multiarm","t12",c11-c10-c01)
   }

   ## (4) random-forest 4-arm T-learner
   arms<-list(c(0,0),c(1,0),c(0,1),c(1,1)); preds<-list()
   for(a in seq_along(arms)){ z<-arms[[a]]; idx<-which(Z1==z[1]&Z2==z[2])
     rf<-randomForest(x=X[idx,,drop=FALSE],y=y[idx],ntree=300); preds[[a]]<-predict(rf,X) }
   av<-anova4(preds[[1]],preds[[2]],preds[[3]],preds[[4]])
   add("RF-4arm-Tlearner","t1",av$t1); add("RF-4arm-Tlearner","t2",av$t2); add("RF-4arm-Tlearner","t12",av$t12)

   ## (5) regularized linear factorial (glmnet S-learner)
   MM<-gdesign(X,Z1,Z2); cvg<-cv.glmnet(MM,y,alpha=0)
   p00<-predict(cvg,gdesign(X,rep(0,nrow(X)),rep(0,nrow(X))),s="lambda.min")[,1]
   p10<-predict(cvg,gdesign(X,rep(1,nrow(X)),rep(0,nrow(X))),s="lambda.min")[,1]
   p01<-predict(cvg,gdesign(X,rep(0,nrow(X)),rep(1,nrow(X))),s="lambda.min")[,1]
   p11<-predict(cvg,gdesign(X,rep(1,nrow(X)),rep(1,nrow(X))),s="lambda.min")[,1]
   gv<-anova4(p00,p10,p01,p11)
   add("Linear-factorial","t1",gv$t1); add("Linear-factorial","t2",gv$t2); add("Linear-factorial","t12",gv$t12)

   do.call(rbind,rows)
  }, error=function(e) data.frame(rep=r,estimator=NA,component=NA,rmse=NA,corr=NA,ate_bias=NA,err=conditionMessage(e)))
  write.csv(out,pf,row.names=FALSE); invisible()
}

R<-15
ncore<-max(1, min(detectCores(), 6))
cat("running",R,"reps on",ncore,"cores\n"); flush.console()
invisible(mclapply(1:R, run_rep, mc.cores=ncore, mc.preschedule=FALSE))

d<-do.call(rbind,lapply(list.files(PARTS,pattern="csv$",full.names=TRUE),function(f)tryCatch(read.csv(f),error=function(e)NULL)))
if("err" %in% names(d)) { errs<-d[!is.na(d$err),]; if(nrow(errs)) {cat("ERRORS in",nrow(errs),"parts:\n"); print(head(unique(errs$err)))} }
d<-d[!is.na(d$estimator),]
agg<-do.call(rbind,lapply(split(d,interaction(d$estimator,d$component,drop=TRUE)),function(s)data.frame(
  estimator=s$estimator[1],component=s$component[1],reps=nrow(s),
  rmse=round(mean(s$rmse),3),rmse_se=round(sd(s$rmse)/sqrt(nrow(s)),3),
  corr=round(mean(s$corr),3),ate_bias=round(mean(s$ate_bias),3))))
saveRDS(list(raw=d,agg=agg),file.path(OUT,"comparison.rds"))
write.csv(agg,file.path(OUT,"comparison_summary.csv"),row.names=FALSE)
cat("\n===== METHOD COMPARISON (RMSE of per-unit effect vs truth) =====\n")
for(cc in c("t1","t2","t12")){cat("\n--",cc,"--\n")
  s<-agg[agg$component==cc,c("estimator","rmse","rmse_se","corr","ate_bias")]
  print(s[order(s$rmse),],row.names=FALSE)}

## figure
if(requireNamespace("ggplot2",quietly=TRUE)){
  library(ggplot2)
  a<-agg; a$component<-factor(a$component,levels=c("t1","t2","t12"),
    labels=c("tau1 (main)","tau2 (main)","tau12 (interaction)"))
  ord<-c("MVBCF-factorial","MVBCF-per-outcome","Separate-forests","grf-multiarm",
         "RF-4arm-Tlearner","Linear-factorial")
  ord<-ord[ord %in% a$estimator]; a$estimator<-factor(a$estimator,levels=rev(ord))
  g<-ggplot(a,aes(rmse,estimator,colour=estimator=="MVBCF-factorial"))+
    geom_errorbarh(aes(xmin=rmse-rmse_se,xmax=rmse+rmse_se),height=.2,colour="grey60")+
    geom_point(size=2.6)+facet_wrap(~component,scales="free_x")+
    scale_colour_manual(values=c("TRUE"="#B2182B","FALSE"="#2166AC"),guide="none")+
    labs(title="Method comparison: RMSE of the per-unit effect vs. truth (lower is better)",
         subtitle=paste0(R," replications; red = proposed factorial MVBCF; grf = real generalized random forest"),
         x="RMSE",y=NULL)+theme_minimal(base_size=9)
  ggsave(file.path(OUT,"comparison_rmse.pdf"),g,width=10,height=3.2)
  ggsave(file.path(OUT,"comparison_rmse.png"),g,width=10,height=3.2,dpi=130)
  cat("figure written\n")
}
cat("\nALL DONE\n")
