suppressMessages({library(Rcpp);library(RcppArmadillo);library(RcppDist);library(ggplot2)})
PKG<-"../mvbcfMT"; sourceCpp(file.path(PKG,"src","mvbcf_multi_engine.cpp"))
for(f in c("components","simulate","fit","predict","bartman","vsup")) source(file.path(PKG,"R",paste0(f,".R")))
set.seed(11)
sim<-simulate_met2(n_gen=24,n_env=6,confounded=TRUE,seed=11); d<-sim$data
NI<-1400
fit<-fit_mvbcf2(d,responses=c("yield_kg_ha","tgw_g","plant_height_cm"),
  treatments=c("Z1","Z2"),gen="gen",env="env",
  n_iter=NI,n_burn=700,keep_every=2,n_tree=80,n_tree_tau=40,
  tree_iters=seq(NI-159,NI),verbose=FALSE)
vars<-fit$covariates
vm_all <- lapply(c("mu","tau_Z1","tau_Z2","tau_Z1:Z2"), function(cl) vivi_matrix(fit, cl, vars))
names(vm_all) <- c("mu","tau_Z1","tau_Z2","tau_Z1:Z2")
saveRDS(list(vm=vm_all, vars=vars), "outputs/vivi_matrices.rds")
mkpanel<-function(cl,ttl) plot_vivi_vsup(fit,cl,vars,title=ttl)
p_mu <-mkpanel("mu","mu (prognostic)")
p_t1 <-mkpanel("tau_Z1","tau_Z1 (drought)")
p_t2 <-mkpanel("tau_Z2","tau_Z2 (second treatment)")
p_t12<-mkpanel("tau_Z1:Z2","tau_Z1:Z2 (interaction)")
leg<-vsup_legend_seq()
suppressMessages(ok<-require(patchwork))
if(ok){ g<-(p_t1|p_t2)/(p_t12|p_mu); ggsave("outputs/vivi_vsup_panel.pdf",g,width=10,height=9)
        ggsave("outputs/vivi_vsup_panel.png",g,width=10,height=9,dpi=110)
        ggsave("outputs/vivi_vsup_legend.pdf",leg,width=3,height=2.6) } else {
        ggsave("outputs/vivi_vsup_tau1.pdf",p_t1,width=5,height=5) }
cat("VIVI panels written\n")
# print the tau1 VIVI value+uncertainty for sanity
vm<-vivi_matrix(fit,"tau_Z1",vars)
cat("tau_Z1 importance (diag):\n"); print(round(sort(diag(vm$value),decreasing=TRUE)[1:5],3))
