library(ggplot2)
library(data.table)

load("FPOP.models.RData")
load("PDPA.models.RData")

FPOP.models[min.log.lambda != -Inf, min(min.log.lambda)]
FPOP.models[max.log.lambda != Inf, max(max.log.lambda)]
## => a grid of values between 5 and 16
exp(seq(5, 16, l=12))

diff.dt <- FPOP.models[, list(
  chunk.name, sample.id,
  fpop.min=min.log.lambda,
  fpop.max=max.log.lambda,
  best.peaks.max,
  n.models,
  fpop.errors=best.errors)][PDPA.targets.dt[, list(
    chunk.name, sample.id,
    pdpa.min=min.log.lambda,
    pdpa.max=max.log.lambda,
    pdpa.errors=errors)], on=list(chunk.name, sample.id)]
stopifnot(nrow(FPOP.models)==nrow(diff.dt))
d <- function(x,y)ifelse(x==y, 0, abs(x-y))
diff.dt[1e-5 < d(pdpa.min, fpop.min) | 1e-5 < d(pdpa.max, fpop.max)]

diff.dt[1e-5 < d(pdpa.max, fpop.max)]
diff.dt[1e-5 < d(pdpa.min, fpop.min)]

diff.dt[pdpa.errors < fpop.errors]
diff.dt[fpop.errors < pdpa.errors]

fpop.lower <- diff.dt[pdpa.min==-Inf & fpop.min!=-Inf]
table(fpop.lower$best.peaks.max) ## from 3 to 119 peaks.
fpop.lower[best.peaks.max==3]
PDPA.modelSelection.dt[chunk.name=="H3K4me3_TDH_immune/15" & sample.id=="McGill0095"]
fpop.lower[best.peaks.max==119]
PDPA.modelSelection.dt[chunk.name=="H3K4me3_PGP_immune/13" & sample.id=="McGill0322"]

diff.dt[, models.diff := 19-n.models]
diff.dt[, errors.diff := pdpa.errors-fpop.errors]
ggplot()+
  geom_point(aes(
    models.diff, errors.diff),
             data=diff.dt,
             shape=1)

ggplot()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(fpop.min, pdpa.min), data=diff.dt, shape=1)

ggplot()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(fpop.max, pdpa.max), data=diff.dt, shape=1)

ggplot()+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_segment(aes(fpop.max, pdpa.max, xend=fpop.min, yend=pdpa.min), data=diff.dt)

both.dt <- diff.dt[, rbind(
  data.table(limit="upper limit", FPOP=fpop.max, PDPA=pdpa.max),
  data.table(limit="lower limit", FPOP=fpop.min, PDPA=pdpa.min))]
gg <- ggplot()+
  ggtitle("PDPA target interval limits often equal FPOP,\nbut PDPA has more infinite limits")+
  theme_bw()+
  facet_grid(. ~ limit)+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_point(aes(FPOP, PDPA), data=both.dt, shape=1)+
  coord_equal()
png("figure-FPOP-PDPA-compare-limits.png", 5, 3, res=100, units="in")
print(gg)
dev.off()

diff.counts <- diff.dt[, list(count=.N), by=list(models.diff, errors.diff)]
more.dt <- data.table(
  x=0,
  hjust=c(1, 0),
  y=2.4,
  label=c("more FPOP ", " more PDPA"))
more.y.dt <- data.table(
  x=-20,
  vjust=c(1.2, -0.2),
  y=0,
  label=c("more FPOP", "more PDPA"))
vline.dt <- data.table(
  xint=c(mean(diff.dt$models.diff), 0))
vline.dt[, label := c(sprintf("mean = %.1f more PDPA models", xint[1]), "")]
gg <- ggplot()+
  ggtitle("FPOP target intervals faster to compute on average, and always more accurate than PDPA")+
  theme_bw()+
  geom_tile(aes(
    models.diff, errors.diff, fill=count),
            data=diff.counts)+
  geom_vline(aes(xintercept=xint), data=vline.dt, color="grey")+
  geom_hline(yintercept=0, color="grey")+
  geom_text(aes(xint, 2, label=paste0(" ", label)), hjust=0, data=vline.dt, color="grey")+
  geom_text(aes(
    models.diff, errors.diff, label=count),
            data=diff.counts)+
  scale_fill_gradient(low="white", high="red")+
  xlab("PDPA models computed - FPOP models computed")+
  ylab("PDPA min errors - \nFPOP min errors")+
  geom_text(aes(x, y, label=label, hjust=hjust), data=more.dt, color="grey")+
  geom_text(aes(x, y, label=label, vjust=vjust), data=more.y.dt, color="grey")
png("figure-FPOP-PDPA-compare.png", 17.5, 2, res=100, units="in")
print(gg)
dev.off()

