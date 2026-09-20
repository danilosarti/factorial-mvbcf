suppressMessages({library(ggplot2)})
r<-readRDS("outputs/robustness.rds"); d<-r$raw
# Panel A: n-sweep convergence (wellspec)
ws<-d[d$scenario=="wellspec",]
A<-do.call(rbind,lapply(split(ws,interaction(ws$n,ws$component,drop=TRUE)),function(s)data.frame(
 n=s$n[1],comp=s$component[1],abs_bias=abs(mean(s$ate_bias)),rmse=mean(s$rmse),
 ate_cov=mean(s$ate_covered),cate_cov=mean(s$cate_cov))))
gA<-ggplot(A,aes(n,rmse,colour=comp))+geom_line()+geom_point()+
  scale_x_log10(breaks=c(300,700,1500,3000))+
  labs(title="A. n-sweep: RMSE -> 0 (consistency)",x="n (log)",y="RMSE of effect",colour=NULL)+
  theme_minimal(base_size=9)+theme(legend.position="bottom")
gAc<-ggplot(A,aes(n,ate_cov,colour=comp))+geom_hline(yintercept=.95,linetype=2,colour="grey55")+
  geom_line()+geom_point()+scale_x_log10(breaks=c(300,700,1500,3000))+ylim(.4,1.02)+
  labs(title="A'. ATE-interval coverage vs n",subtitle="dashed=0.95; noisy (12 reps) but ~nominal at large n",x="n (log)",y="coverage",colour=NULL)+
  theme_minimal(base_size=9)+theme(legend.position="none")
# Panel B: coverage by scenario (n=1200/1500 main effects), CATE vs ATE
sc<-d[d$n %in% c(1200) & d$component %in% c("tau_Z1","tau_Z2"),]
B<-do.call(rbind,lapply(split(sc,interaction(sc$scenario,sc$component,drop=TRUE)),function(s)data.frame(
 scenario=s$scenario[1],comp=s$component[1],
 CATE=mean(s$cate_cov),ATE=mean(s$ate_covered))))
Bl<-reshape(B,varying=c("CATE","ATE"),v.names="cov",timevar="type",times=c("CATE","ATE"),direction="long")
Bl$scenario<-factor(Bl$scenario,levels=c("strong_confounding","heavytail","T3_correct","T3_underspec","unmeasured"))
gB<-ggplot(Bl,aes(cov,scenario,colour=type,shape=comp))+
  geom_vline(xintercept=.95,linetype=2,colour="grey55")+
  geom_point(size=2.6,position=position_dodge(width=.5))+xlim(0,1.02)+
  scale_colour_manual(values=c(ATE="#B2182B",CATE="#2166AC"))+
  labs(title="B. Coverage under misspecification (main effects, n=1200)",
       subtitle="unmeasured confounding collapses to ~0 (correctly); others near/below nominal",
       x="empirical 95% coverage",y=NULL,colour=NULL,shape=NULL)+
  theme_minimal(base_size=9)+theme(legend.position="bottom")
suppressMessages({ok<-require(patchwork)})
if(ok){ g<-(gA|gAc)/gB+plot_layout(heights=c(1,1.1)) } else g<-gB
ggsave("outputs/robustness_summary.pdf",g,width=9,height=7)
ggsave("outputs/robustness_summary.png",g,width=9,height=7,dpi=120)
cat("figure written\n")
