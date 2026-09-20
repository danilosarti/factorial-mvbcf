suppressMessages(library(ggplot2))
r <- readRDS("outputs/mc_validation.rds"); a <- r$agg
a$component <- factor(a$component, levels=c("tau_Z1","tau_Z2","tau_Z1:Z2"))
# panel 1: relative bias +/- 2 MC-SE
a$rb <- a$rel_bias_pct; a$rse <- 100*a$mc_se_bias/ a$mean_rmse # approx scale; recompute properly:
# use ate_bias mc-se relative to |true ate|; approximate via raw
raw <- r$raw
bias <- do.call(rbind, lapply(split(raw, raw$component), function(s) data.frame(
  component=s$component[1],
  rb=100*mean(s$ate_bias)/mean(abs(s$ate_true)),
  rse=100*(sd(s$ate_bias)/sqrt(nrow(s)))/mean(abs(s$ate_true)))))
bias$component <- factor(bias$component, levels=c("tau_Z1","tau_Z2","tau_Z1:Z2"))
g1 <- ggplot(bias, aes(rb, component)) +
  geom_vline(xintercept=0, linetype=2, colour="grey55") +
  geom_errorbarh(aes(xmin=rb-2*rse, xmax=rb+2*rse), height=.15, colour="grey45") +
  geom_point(colour="#2166AC", size=2.6) +
  labs(title="ATE relative bias (%) +/- 2 MC-SE", subtitle="straddles 0 = unbiased", x="relative bias (%)", y=NULL) +
  theme_minimal(base_size=10)
# panel 2: coverage
cov <- rbind(data.frame(component=a$component, type="CATE (per-unit)", cov=a$cate_coverage),
             data.frame(component=a$component, type="ATE interval", cov=a$ate_interval_coverage))
g2 <- ggplot(cov, aes(cov, component, colour=type)) +
  geom_vline(xintercept=0.95, linetype=2, colour="grey55") +
  geom_point(size=2.6, position=position_dodge(width=.4)) +
  scale_colour_manual(values=c("ATE interval"="#B2182B","CATE (per-unit)"="#2166AC")) +
  labs(title="95% credible-interval coverage", subtitle="dashed = nominal 0.95", x="empirical coverage", y=NULL, colour=NULL) +
  xlim(0.7,1) + theme_minimal(base_size=10) + theme(legend.position="bottom")
suppressMessages(library(patchwork))
if("patchwork" %in% rownames(installed.packages())){ g <- g1+g2 } else { g <- g1 }
ggsave("outputs/mc_validation_bias_coverage.pdf", g, width=9, height=3)
ggsave("outputs/mc_validation_bias_coverage.png", g, width=9, height=3, dpi=130)
cat("figure written\n")
