suppressMessages(library(ggplot2))
r<-readRDS("outputs/comparison.rds"); a<-r$agg
a$component<-factor(a$component,levels=c("t1","t2","t12"),
  labels=c("tau1 (main)","tau2 (main)","tau12 (interaction)"))
ord<-c("MVBCF-factorial","MVBCF-per-outcome","Separate-forests","RF-4arm-Tlearner","Linear-factorial")
a$estimator<-factor(a$estimator,levels=rev(ord))
g<-ggplot(a,aes(rmse,estimator,colour=estimator=="MVBCF-factorial"))+
  geom_errorbarh(aes(xmin=rmse-rmse_se,xmax=rmse+rmse_se),height=.2,colour="grey60")+
  geom_point(size=2.6)+facet_wrap(~component,scales="free_x")+
  scale_colour_manual(values=c("TRUE"="#B2182B","FALSE"="#2166AC"),guide="none")+
  labs(title="Method comparison: RMSE of the per-unit effect vs. truth (lower is better)",
       subtitle="15 replications; red = proposed factorial MVBCF. Separate-forests cannot estimate the interaction.",
       x="RMSE",y=NULL)+theme_minimal(base_size=9)
ggsave("outputs/comparison_rmse.pdf",g,width=10,height=3.2)
ggsave("outputs/comparison_rmse.png",g,width=10,height=3.2,dpi=130)
cat("fig written\n")
