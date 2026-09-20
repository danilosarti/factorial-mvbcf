# Cheap re-plot of the 3-treatment vivid VIVI from saved matrices (no refit).
suppressMessages({library(ggplot2);library(patchwork)})
PKG<-"../mvbcfMT"; source(file.path(PKG,"R","vsup.R"))
R<-readRDS("outputs/vivi3_matrices.rds")
hm_from<-function(cl,vars,ttl,cvmax=2,max_desat=0.85){
  vm<-R$vm[[cl]]; V<-length(vars)
  df<-expand.grid(ri=1:V,ci=1:V)
  df$value<-mapply(function(r,c)vm$value[r,c],df$ri,df$ci)
  df$sd<-mapply(function(r,c)vm$uncertainty[r,c],df$ri,df$ci)
  df$diag<-df$ri==df$ci; df$cv<-pmin(df$sd/pmax(df$value,1e-6),cvmax)/cvmax; df$w<-df$cv*max_desat
  vmd<-max(df$value[df$diag]); vmo<-max(df$value[!df$diag])
  di<-df$diag; oi<-!df$diag; df$hex<-NA
  df$hex[di]<-mvbcfMT_vsupmix(df$value[di]/max(vmd,1e-9),df$w[di],"#FFFFCC","#08306B","#BFB8AE")
  df$hex[oi]<-mvbcfMT_vsupmix(df$value[oi]/max(vmo,1e-9),df$w[oi],"#FFFFCC","#67000D","#BFB8AE")
  df$row<-factor(vars[df$ri],levels=rev(vars)); df$col<-factor(vars[df$ci],levels=vars)
  hm<-ggplot(df,aes(col,row,fill=hex))+geom_tile(colour="grey88",linewidth=.3)+scale_fill_identity()+
    coord_equal()+labs(title=ttl,x=NULL,y=NULL)+theme_minimal(base_size=8)+
    theme(axis.text.x=element_text(angle=45,hjust=1),axis.text.y=element_text(size=6),panel.grid=element_blank())
  attr(hm,"maxes")<-c(imp=vmd,int=vmo); hm }
mvbcfMT_vsupmix<-get(".vsup_mix")
cl<-"mu"; vars<-c(R$covariates,paste0("pro_",R$treatments))
hm<-hm_from(cl,vars,"VIVI (VSUP): mu -- 3 treatments"); mx<-attr(hm,"maxes")
g<-hm | (vsup_wedge(mx["int"],"int")/vsup_wedge(mx["imp"],"imp"))
ggsave("outputs/vivi3_vivid_mu.pdf",g,width=9.5,height=5.2); ggsave("outputs/vivi3_vivid_mu.png",g,width=9.5,height=5.2,dpi=120)
cat("replotted\n")
