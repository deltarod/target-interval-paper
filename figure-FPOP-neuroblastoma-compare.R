source("packages.R")

load("FPOP.neuroblastoma.RData")

diff.dt <- dcast(
  FPOP.neuroblastoma,
  profile.id + chromosome ~ algo,
  value.var=c("min.log.lambda", "max.log.lambda", "errors", "n.models"))

d <- function(x,y)ifelse(x==y, 0, abs(x-y))
diff.dt[1e-5 < d(max.log.lambda_PDPA, max.log.lambda_proposed)]
diff.dt[1e-5 < d(min.log.lambda_PDPA, min.log.lambda_proposed)]

## Same number of errors in all methods.
FPOP.neuroblastoma[errors != 0]
diff.dt[errors_proposed < errors_grid.search] 
diff.dt[errors_grid.search < errors_proposed] 
diff.dt[errors_proposed < errors_PDPA] 
diff.dt[errors_PDPA < errors_proposed]

diff.dt[, table(n.models_proposed)]

taller.dt <- melt(FPOP.neuroblastoma, measure.vars=c("min.log.lambda", "max.log.lambda"))
both.dt <- dcast(taller.dt, profile.id + chromosome + variable ~ algo)
gg <- ggplot()+
  theme_bw()+
  facet_grid(. ~ variable)+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_abline(slope=1, intercept=0, color="grey")+
  ##geom_vline(aes(xintercept=log.penalty), color="red", data=data.table(log.penalty=log(as.numeric(penalty.vec))))+
  scale_x_continuous(breaks=log(as.numeric(penalty.vec)))+
  scale_y_continuous(breaks=log(as.numeric(penalty.vec)))+
  geom_point(aes(proposed, grid.search), data=both.dt, shape=1)+
  coord_equal()
png("figure-FPOP-gridsearch-compare-limits.png", 9, 5, res=100, units="in")
print(gg)
dev.off()

diff.counts <- diff.dt[, list(count=.N), by=list(models.diff, errors.diff)]
more.dt <- data.table(
  x=c(-1, 1),
  hjust=c(1, 0),
  y=2.5,
  label=c("<- more proposed", "more grid search ->"))
more.y.dt <- data.table(
  x=-25,
  vjust=c(1.2, -0.3),
  y=0,
  label=c("more proposed", "more grid search"))
vline.dt <- data.table(
  xint=c(mean(diff.dt$models.diff, na.rm=TRUE), 0))
vline.dt[, label := c(sprintf("mean = %.1f more grid search models", xint[1]), "")]
gg <- ggplot()+
  ggtitle("Proposed target interval search faster to compute on average, and always more accurate than grid search")+
  theme_bw()+
  geom_tile(aes(
    models.diff, errors.diff, fill=count),
            data=diff.counts)+
  geom_vline(aes(xintercept=xint), data=vline.dt, color="grey")+
  geom_hline(yintercept=0, color="grey")+
  geom_text(aes(xint, 3.5, label=paste0(" ", label)), hjust=0, data=vline.dt, color="grey")+
  geom_text(aes(
    models.diff, errors.diff, label=count),
            data=diff.counts)+
  scale_fill_gradient(low="white", high="red")+
  xlab("Difference in models computed (Grid search - Proposed target interval search)")+
  ylab("Difference in min errors\n(Grid search - Proposed target interval search)")+
  geom_text(aes(x, y, label=label, hjust=hjust), data=more.dt, color="grey")+
  geom_text(aes(x, y, label=label, vjust=vjust), data=more.y.dt, color="grey")
png("figure-FPOP-gridsearch-compare.png", 17.5, 4, res=100, units="in")
print(gg)
dev.off()

