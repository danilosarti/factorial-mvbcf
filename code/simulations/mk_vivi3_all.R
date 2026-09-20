suppressMessages({library(ggplot2);library(patchwork)})
PKG<-"../mvbcfMT"; source(file.path(PKG,"R","vsup.R"))
mix<-get(".vsup_mix")
R<-readRDS("outputs/vivi3_matrices.rds")
hm<-function(cl,ttl,cvmax=2,max_desat=0.85,ytext=6){
  vm<-R$vm[[cl]]; vars<-rownames(vm$value); V<-length(vars)
  df<-expand.grid(ri=1:V,ci=1:V)
  df$value<-mapply(function(r,c)vm$value[r,c],df$ri,df$ci)
  df$sd<-mapply(function(r,c)vm$uncertainty[r,c],df$ri,df$ci)
  df$diag<-df$ri==df$ci; df$cv<-pmin(df$sd/pmax(df$value,1e-6),cvmax)/cvmax; df$w<-df$cv*max_desat
  vmd<-max(df$value[df$diag]); vmo<-max(df$value[!df$diag])
  di<-df$diag; oi<-!df$diag; df$hex<-NA
  df$hex[di]<-mix(df$value[di]/max(vmd,1e-9),df$w[di],"#FFFFCC","#08306B","#BFB8AE")
  df$hex[oi]<-mix(df$value[oi]/max(vmo,1e-9),df$w[oi],"#FFFFCC","#67000D","#BFB8AE")
  df$row<-factor(vars[df$ri],levels=rev(vars)); df$col<-factor(vars[df$ci],levels=vars)
  g<-ggplot(df,aes(col,row,fill=hex))+geom_tile(colour="grey88",linewidth=.25)+scale_fill_identity()+
    coord_equal()+labs(title=ttl,subtitle=sprintf("max Vimp=%.3f  max Vint=%.3f",vmd,vmo),x=NULL,y=NULL)+
    theme_minimal(base_size=7)+theme(axis.text.x=element_text(angle=45,hjust=1,size=5),
      axis.text.y=element_text(size=ytext),panel.grid=element_blank(),
      plot.subtitle=element_text(size=6,colour="grey40"))
  attr(g,"maxes")<-c(imp=vmd,int=vmo); g }

# individual full figures (heatmap + own wedges) for all 6 effect/interaction forests
comps<-c("tau_Z1","tau_Z2","tau_Z3","tau_Z1:Z2","tau_Z1:Z3","tau_Z2:Z3")
titles<-c("tau_Z1 (treatment 1)","tau_Z2 (treatment 2)","tau_Z3 (treatment 3)",
          "tau_Z1:Z2 (interaction)","tau_Z1:Z3 (interaction)","tau_Z2:Z3 (interaction)")
for(i in seq_along(comps)){ h<-hm(comps[i],titles[i],ytext=7); mx<-attr(h,"maxes")
  fig<- h | (vsup_wedge(mx["int"],"int")/vsup_wedge(mx["imp"],"imp"))
  fn<-gsub("[:]","",comps[i]); ggsave(sprintf("outputs/vivi3_vivid_%s.pdf",fn),fig,width=8.5,height=4.6)
  ggsave(sprintf("outputs/vivi3_vivid_%s.png",fn),fig,width=8.5,height=4.6,dpi=110) }

# gallery: 6 heatmaps (2x3) + one reference fan-legend pair (encoding; scales per forest)
hs<-lapply(seq_along(comps),function(i)hm(comps[i],titles[i]))
leg<-(vsup_wedge(1,"int",title="Vint")/vsup_wedge(1,"imp",title="Vimp"))
gallery<-(wrap_plots(hs,ncol=3)) | leg
gallery<-gallery + plot_layout(widths=c(6,1))
ggsave("outputs/vivi3_gallery.pdf",gallery,width=15,height=8)
ggsave("outputs/vivi3_gallery.png",gallery,width=15,height=8,dpi=95)
cat("all vivid figures written\n")
