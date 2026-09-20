suppressMessages({library(ggplot2)})
PKG<-"../mvbcfMT"; source(file.path(PKG,"R","vsup.R"))  # palette helpers (no engine needed)
R<-readRDS("outputs/vivi_matrices.rds"); vars<-R$vars
panel<-function(cl,ttl){ vm<-R$vm[[cl]]
  df<-expand.grid(row=vars,col=vars,stringsAsFactors=FALSE)
  df$value<-mapply(function(r,c)vm$value[r,c],df$row,df$col)
  df$unc<-mapply(function(r,c)vm$uncertainty[r,c],df$row,df$col)
  df$hex<-vsup_palette_seq(df$value,df$unc)
  df$row<-factor(df$row,levels=rev(vars)); df$col<-factor(df$col,levels=vars)
  ggplot(df,aes(col,row,fill=hex))+geom_tile(colour="white",linewidth=.3)+
    scale_fill_identity()+coord_equal()+
    labs(title=ttl,subtitle="diag=importance, off-diag=interaction; grey=uncertain",x=NULL,y=NULL)+
    theme_minimal(base_size=8)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())}
p1<-panel("tau_Z1","tau_Z1 (drought)"); p2<-panel("tau_Z2","tau_Z2 (second treatment)")
p3<-panel("tau_Z1:Z2","tau_Z1:Z2 (interaction)"); p4<-panel("mu","mu (prognostic)")
leg<-vsup_legend_seq()
if(require(patchwork,quietly=TRUE)){ g<-(p1|p2)/(p3|p4)
  ggsave("outputs/vivi_vsup_panel.pdf",g,width=10,height=9)
  ggsave("outputs/vivi_vsup_panel.png",g,width=10,height=9,dpi=110)
  ggsave("outputs/vivi_vsup_legend.pdf",leg,width=3,height=2.6)
  ggsave("outputs/vivi_vsup_legend.png",leg,width=3,height=2.6,dpi=120)}
cat("replotted\n")
