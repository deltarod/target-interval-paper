library(data.table)
library(parallel)
library(PeakSegPipeline)
counts.RData.vec <- Sys.glob("data/H*/*/counts.RData")
i.vec <- seq_along(counts.RData.vec)
##i.vec <- 1:2
penalty.vec <- paste(exp(seq(5, 16, l=12)))

FPOP.gridsearch.list <- list()
for(counts.RData.i in i.vec){
  counts.RData <- counts.RData.vec[[counts.RData.i]]
  chunk.dir <- dirname(counts.RData)
  chunk.name <- sub("data/", "", chunk.dir)
  set.name <- dirname(chunk.name)
  chunk.id <- basename(chunk.name)
  PeakSegFPOP_gridsearch.rds <- file.path(chunk.dir, "PeakSegFPOP_gridsearch.rds")
  if(file.exists(PeakSegFPOP_gridsearch.rds)){
    target.dt <- readRDS(PeakSegFPOP_gridsearch.rds)
  }else{
    cat(sprintf(
      "%4d / %4d %s\n",
      counts.RData.i,
      length(counts.RData.vec),
      chunk.dir))
    load(counts.RData)
    regions.RData <- file.path(chunk.dir, "regions.RData")
    load(regions.RData)
    regions.dt <- data.table(regions)
    setkey(regions.dt, cell.type, sample.id)
    counts.dt <- data.table(counts)
    target.dt <- counts.dt[, {
      chrom <- paste(regions.dt$chrom[1])
      m <- data.table(sample.id)
      sample.regions <- regions.dt[m, on=list(sample.id)]
      problem <- data.table(
        chrom,
        problemStart=min(chromStart),
        problemEnd=max(chromEnd))
      sample.dir <- file.path(
        chunk.dir, "PeakSegFPOP_samples", sample.id,
        problem[, sprintf("%s:%d-%d", chrom, problemStart, problemEnd)])
      loss.dt <- do.call(rbind, lapply(penalty.vec, function(penalty){
        fit <- problem.PeakSegFPOP(sample.dir, penalty)
        fit$error.regions <- PeakError::PeakErrorChrom(fit$segments[status=="peak"], sample.regions)
        fit$loss[, errors := with(fit$error.regions, sum(fn+fp)) ]
        with(fit, loss[timing, on=list(penalty)])
      }))
      ms.dt <- data.table(penaltyLearning::modelSelection(loss.dt, "total.cost", "peaks"))
      target <- penaltyLearning::targetIntervals(ms.dt, "bases")
      min.err.dt <- ms.dt[errors==min(errors)]
      data.table(
        n.models=nrow(loss.dt),
        best.errors=target$errors,
        best.peaks.min=min(min.err.dt$peaks),
        best.peaks.max=max(min.err.dt$peaks),
        min.log.lambda=target$min.log.lambda,
        max.log.lambda=target$max.log.lambda)
    }, by=list(cell.type, sample.id)]
    saveRDS(target.dt, PeakSegFPOP_gridsearch.rds)
  }
  data.table(
    chunk.name, chunk.id, set.name,
    target.dt)
}
FPOP.gridsearch <- do.call(rbind, FPOP.gridsearch.list)
 
save(FPOP.gridsearch, file="FPOP.gridsearch.RData")
