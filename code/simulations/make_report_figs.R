# Master analysis: generates every figure/number used in the report.
suppressMessages({library(Rcpp); library(RcppArmadillo); library(RcppDist); library(ggplot2)})
PKG <- c("../mvbcfMT","mvbcfMT")[dir.exists(c("../mvbcfMT","mvbcfMT"))][1]
OUT <- "outputs"; dir.create(OUT, showWarnings = FALSE)
sourceCpp(file.path(PKG,"src","mvbcf_multi_engine.cpp"))
for (f in c("components","simulate","fit","predict","bartman","vsup"))
  source(file.path(PKG,"R",paste0(f,".R")))
cr <- function(a,b) cor(as.vector(a), as.vector(b))
res <- list()

## ============ PART A: simulated 2-treatment MET (truth known) ==============
cat("[A] fitting 2-treatment MET ...\n")
set.seed(11)
sim <- simulate_met2(n_gen = 24, n_env = 6, confounded = TRUE, seed = 11)
d <- sim$data
fitA <- fit_mvbcf2(d, responses = c("yield_kg_ha","tgw_g","plant_height_cm"),
                   treatments = c("Z1","Z2"), gen = "gen", env = "env",
                   n_iter = 1200, n_burn = 600, keep_every = 2,
                   n_tree = 80, n_tree_tau = 40, verbose = FALSE)
tru <- function(p) as.matrix(d[, paste0("true_",p,c("_yield","_tgw","_height"))])
est <- function(l) component_effect(fitA, l, test = FALSE)$mean
recov <- data.frame(
  comp = c("tau1 (drought)","tau2 (inoculant)","tau12 (synergy)"),
  corr = c(cr(tru("tau1"),est("tau_Z1")), cr(tru("tau2"),est("tau_Z2")), cr(tru("tau12"),est("tau_Z1:Z2"))),
  ate_true = c(mean(tru("tau1")[,1]),mean(tru("tau2")[,1]),mean(tru("tau12")[,1])),
  ate_est  = c(mean(est("tau_Z1")[,1]),mean(est("tau_Z2")[,1]),mean(est("tau_Z1:Z2")[,1])))
res$recovA <- recov; print(recov)
# coverage
cvg <- sapply(c("tau_Z1","tau_Z2","tau_Z1:Z2"), function(l){
  ce <- component_effect(fitA,l,test=FALSE); tp <- switch(l,"tau_Z1"="tau1","tau_Z2"="tau2","tau_Z1:Z2"="tau12")
  mean(tru(tp) >= ce$lo & tru(tp) <= ce$hi)})
res$coverage <- cvg; cat("coverage(all traits):", paste(names(cvg),round(cvg,2)), "\n")

# --- Fig A1: recovery scatter ---
pd <- rbind(data.frame(component="tau1 (drought)", true=tru("tau1")[,1], est=est("tau_Z1")[,1]),
            data.frame(component="tau2 (inoculant)", true=tru("tau2")[,1], est=est("tau_Z2")[,1]),
            data.frame(component="tau12 (synergy)", true=tru("tau12")[,1], est=est("tau_Z1:Z2")[,1]))
lab <- do.call(rbind, lapply(split(pd,pd$component), function(s)
  data.frame(component=s$component[1], r=cr(s$true,s$est),
             x=min(s$true), y=max(s$est))))
gA1 <- ggplot(pd,aes(true,est))+geom_point(alpha=.35,colour="#2166AC",size=.8)+
  geom_abline(slope=1,intercept=0,linetype=2,colour="grey35")+
  geom_text(data=lab,aes(x,y,label=sprintf("r=%.2f",r)),hjust=0,vjust=1,size=3)+
  facet_wrap(~component,scales="free",nrow=1)+
  labs(title="A. Two-treatment MVBCF recovers the known effect surfaces (yield)",
       x="true effect",y="posterior-mean effect")+theme_minimal(base_size=9)
ggsave(file.path(OUT,"figA1_recovery.pdf"),gA1,width=9,height=3)

# --- Fig A2: per-component importance (correct drivers) ---
gA2 <- plot_component_importance(fitA, top=6)
ggsave(file.path(OUT,"figA2_importance.pdf"),gA2,width=8.5,height=6)
res$vimp_tau1 <- head(vimp_table(fitA,"tau_Z1"),5)
res$vimp_tau2 <- head(vimp_table(fitA,"tau_Z2"),5)
res$vimp_tau12<- head(vimp_table(fitA,"tau_Z1:Z2"),5)
cat("drivers tau1:", paste(res$vimp_tau1$variable,collapse=", "),"\n")
cat("drivers tau2:", paste(res$vimp_tau2$variable,collapse=", "),"\n")

# --- Fig A3: VSUP maps of tau1 (drought) and tau12 (interaction) ---
v1  <- plot_effect_vsup(fitA,"tau_Z1","yield_kg_ha")
v12 <- plot_effect_vsup(fitA,"tau_Z1:Z2","yield_kg_ha")
leg <- vsup_legend()
ggsave(file.path(OUT,"figA3a_vsup_tau1.pdf"),v1,width=5,height=5.2)
ggsave(file.path(OUT,"figA3b_vsup_tau12.pdf"),v12,width=5,height=5.2)
ggsave(file.path(OUT,"figA3c_vsup_legend.pdf"),leg,width=3.4,height=3)
cat("[A] figures done\n")

## ============ PART B: REAL data - npk factorial ============================
cat("[B] fitting real npk factorial ...\n")
data(npk)
npk2 <- npk
npk2$N <- as.integer(npk2$N==1); npk2$P <- as.integer(npk2$P==1); npk2$K <- as.integer(npk2$K==1)
# block one-hot as covariates
for (b in levels(npk$block)) npk2[[paste0("blk_",b)]] <- as.integer(npk$block==b)
covs <- grep("^blk_",names(npk2),value=TRUE)
fitB <- fit_mvbcf_multi(npk2, responses="yield", treatments=c("N","P","K"),
                        covariates=covs, order=2,
                        n_iter=3000, n_burn=1500, keep_every=3,
                        n_tree=40, n_tree_tau=20, min_nodesize=2, verbose=FALSE)
# main/interaction effects via component_effect (each tau_S surface; ATE = mean)
mains <- lapply(c("tau_N","tau_P","tau_K","tau_N:P","tau_N:K","tau_P:K"), function(l){
  ce<-component_effect(fitB,l,test=FALSE); c(l, mean(ce$mean), mean(ce$lo), mean(ce$hi))})
bt <- as.data.frame(do.call(rbind, mains)); names(bt)<-c("effect","ate","lo","hi")
bt[,2:4] <- lapply(bt[,2:4], as.numeric)
res$npk <- bt; cat("\n--- npk MVBCF effects (yield) ---\n"); print(bt, digits=3, row.names=FALSE)
# classical benchmark
av <- summary(aov(yield ~ block + N*P*K, data=npk))
res$npk_aov <- av; cat("\nClassical aov(yield~block+N*P*K):\n"); print(av)
# lm coefficients for N,P,K (main) for direct comparison
lmfit <- lm(yield ~ N+P+K, data=npk2); res$npk_lm <- coef(summary(lmfit))
cat("\nlm main-effect coefs:\n"); print(round(coef(summary(lmfit)),2))
# Fig B: effect estimates with CI
bt$effect <- factor(bt$effect, levels=rev(bt$effect))
gB <- ggplot(bt, aes(ate, effect))+
  geom_vline(xintercept=0,colour="grey60",linetype=2)+
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=.2,colour="grey50")+
  geom_point(colour="#B2182B",size=2.4)+
  labs(title="B. Real data (npk): MVBCF effects on yield with 95% CrI",
       subtitle="main effects N,P,K and pairwise interactions; dashed = 0",
       x="average effect on yield",y=NULL)+theme_minimal(base_size=10)
ggsave(file.path(OUT,"figB_npk_effects.pdf"),gB,width=7,height=3.4)
cat("[B] done\n")

## ============ PART C: T=3 recovery incl. spurious interaction ==============
cat("[C] T=3 recovery figure ...\n")
set.seed(7)
simC <- simulate_multi(n=1200, p=6, q=2, Ti=3, seed=7)
fitC <- fit_mvbcf_multi(simC$data, responses=simC$responses, treatments=simC$treatments,
                        covariates=simC$covariates, order=2,
                        n_iter=700, n_burn=350, keep_every=2, n_tree=40, n_tree_tau=25, verbose=FALSE)
mn <- sapply(1:3, function(t) cr(simC$truth$main[[t]], component_effect(fitC,paste0("tau_Z",t),test=FALSE)$mean))
i12 <- cr(simC$truth$pair[["1:2"]], component_effect(fitC,"tau_Z1:Z2",test=FALSE)$mean)
i23 <- cr(simC$truth$pair[["2:3"]], component_effect(fitC,"tau_Z2:Z3",test=FALSE)$mean)
sd_sp <- sd(component_effect(fitC,"tau_Z1:Z3",test=FALSE)$mean)
sd_re <- mean(c(sd(component_effect(fitC,"tau_Z1:Z2",test=FALSE)$mean),
                sd(component_effect(fitC,"tau_Z2:Z3",test=FALSE)$mean)))
res$T3 <- list(main=mn, i12=i12, i23=i23, sd_spurious=sd_sp, sd_real=sd_re)
cat(sprintf("main corr: %.2f %.2f %.2f | i12=%.2f i23=%.2f | sd spurious/real=%.2f/%.2f\n",
            mn[1],mn[2],mn[3],i12,i23,sd_sp,sd_re))
dfC <- data.frame(
  effect=c("Z1","Z2","Z3","Z1:Z2\n(real)","Z2:Z3\n(real)","Z1:Z3\n(spurious)"),
  metric=c(mn,i12,i23, sd_sp/sd_re),
  kind=c(rep("main",3),"real int","real int","spurious"))
# for the spurious one plot the shrinkage ratio instead of corr; annotate
res$dfC <- dfC
saveRDS(res, file.path(OUT,"report_results.rds"))
cat("[C] done; all results saved.\nALL DONE\n")
