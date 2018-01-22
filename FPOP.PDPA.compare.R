library(data.table)

load("FPOP.models.RData")
load("PDPA.models.RData")

diff.dt <- FPOP.models[, list(
  chunk.name, sample.id,
  fpop.min=min.log.lambda,
  fpop.max=max.log.lambda,
  fpop.errors=best.errors)][PDPA.targets.dt[, list(
    chunk.name, sample.id,
    pdpa.min=min.log.lambda,
    pdpa.max=max.log.lambda,
    pdpa.errors=errors)], on=list(chunk.name, sample.id)]
stopifnot(nrow(FPOP.models)==nrow(diff.dt))
d <- function(x,y)ifelse(x==y, 0, abs(x-y))
diff.dt[1e-5 < d(pdpa.min, fpop.min) | 1e-5 < d(pdpa.max, fpop.max)]
