suppressMessages({library(Rcpp);library(RcppArmadillo);library(RcppDist);library(ggplot2);library(patchwork)})
PKG<-"../mvbcfMT"; sourceCpp(file.path(PKG,"src","mvbcf_multi_engine.cpp"))
for(f in c("components","simulate","fit","predict","bartman","vsup")) source(file.path(PKG,"R",paste0(f,".R")))
set.seed(7)
sim<-simulate_multi(n=1300,p=6,q=2,Ti=3,seed=7)
metnames<-c("temp_mean","precip_total","vpd","radiation","soil_OM","soil_N")
names(sim$data)[match(sim$covariates,names(sim$data))]<-metnames; sim$covariates<-metnames
NI<-1100
fit<-fit_mvbcf_multi(sim$data,responses=sim$responses,treatments=sim$treatments,
  covariates=sim$covariates,order=2,n_iter=NI,n_burn=550,keep_every=2,
  n_tree=50,n_tree_tau=22,tree_iters=seq(NI-149,NI),verbose=FALSE)
cat("components:",paste(fit$components$labels,collapse=", "),"\n")
# save matrices for cheap replot
vm_all<-lapply(fit$components$labels,function(cl){
  vv<-if(cl=="mu") c(fit$covariates, paste0("pro_",fit$treatments)) else fit$covariates
  vivi_matrix(fit,cl,vv)}); names(vm_all)<-fit$components$labels
saveRDS(list(vm=vm_all,covariates=fit$covariates,treatments=fit$treatments,labels=fit$components$labels),
        "outputs/vivi3_matrices.rds")

flagship<-function(cl,ttl){
  vv<-if(cl=="mu") c(fit$covariates, paste0("pro_",fit$treatments)) else fit$covariates
  hm<-plot_vivi_vivid(fit,cl,vv,title=ttl); mx<-attr(hm,"maxes")
  wint<-vsup_wedge(mx["int"],"int"); wimp<-vsup_wedge(mx["imp"],"imp")
  hm | (wint/wimp) + patchwork::plot_layout(heights=c(1,1))
}
g_mu <- flagship("mu","VIVI (VSUP): mu (prognostic) -- 3 treatments")
g_int<- flagship("tau_Z2:Z3","VIVI (VSUP): tau_Z2:Z3 (interaction)")
ggsave("outputs/vivi3_vivid_mu.pdf", g_mu, width=9.5,height=5.2)
ggsave("outputs/vivi3_vivid_mu.png", g_mu, width=9.5,height=5.2,dpi=120)
ggsave("outputs/vivi3_vivid_interaction.pdf", g_int, width=9.5,height=5.2)
ggsave("outputs/vivi3_vivid_interaction.png", g_int, width=9.5,height=5.2,dpi=120)

# 7-component small panel (simple sequential VSUP) to show 3-treatment scaling
panels<-lapply(fit$components$labels,function(cl){
  vv<-if(cl=="mu") c(fit$covariates, paste0("pro_",fit$treatments)) else fit$covariates
  plot_vivi_vsup(fit,cl,vv,title=cl)+theme(plot.subtitle=element_blank())})
G<-wrap_plots(panels,ncol=4)
ggsave("outputs/vivi3_panel_all.pdf",G,width=14,height=7)
ggsave("outputs/vivi3_panel_all.png",G,width=14,height=7,dpi=95)
cat("DONE vivi3\n")
