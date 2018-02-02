source("packages.R")
counts.RData.vec <- Sys.glob("data/H*/*/counts.RData")
i.vec <- seq_along(counts.RData.vec)
##i.vec <- 1:2
FPOP.models.list <- mclapply(i.vec, function(counts.RData.i){
  counts.RData <- counts.RData.vec[[counts.RData.i]]
  chunk.dir <- dirname(counts.RData)
  chunk.name <- sub("data/", "", chunk.dir)
  set.name <- dirname(chunk.name)
  chunk.id <- basename(chunk.name)
  PeakSegFPOP_targets.rds <- file.path(chunk.dir, "PeakSegFPOP_targets.rds")
  if(file.exists(PeakSegFPOP_targets.rds)){
    target.dt <- readRDS(PeakSegFPOP_targets.rds)
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
      problem <- data.table(
        chrom,
        problemStart=min(chromStart),
        problemEnd=max(chromEnd))
      sample.dir <- file.path(
        chunk.dir, "PeakSegFPOP_samples", sample.id,
        problem[, sprintf("%s:%d-%d", chrom, problemStart, problemEnd)])
      dir.create(sample.dir, showWarnings=FALSE, recursive=TRUE)
      prob.cov <- data.table(chrom, chromStart, chromEnd, coverage)
      fwrite(
        prob.cov,
        file.path(sample.dir, "coverage.bedGraph"),
        quote=FALSE,
        sep="\t",
        col.names=FALSE)
      fwrite(
        problem,
        file.path(sample.dir, "problem.bed"),
        quote=FALSE,
        sep="\t",
        col.names=FALSE)
      meta.dt <- data.table(cell.type, sample.id)
      labels <- regions.dt[meta.dt]
      fwrite(
        labels[, list(chrom, chromStart, chromEnd, annotation)],
        file.path(sample.dir, "labels.bed"),
        quote=FALSE,
        sep="\t",
        col.names=FALSE)
      print(sample.dir)
      L <- problem.target(sample.dir)
      print(L$models)
      errors.dt <- L$models[, errors := fp + fn]
      min.err.dt <- errors.dt[min(errors)==errors]
      data.table(
        n.models=nrow(L$models),
        best.errors=min.err.dt$errors[1],
        best.peaks.min=min(min.err.dt$peaks),
        best.peaks.max=max(min.err.dt$peaks),
        min.log.lambda=L$target[1],
        max.log.lambda=L$target[2])
    }, by=list(cell.type, sample.id)]
    saveRDS(target.dt, PeakSegFPOP_targets.rds)
  }
  data.table(
    chunk.name, chunk.id, set.name,
    target.dt)
})
FPOP.models <- do.call(rbind, FPOP.models.list)
 
save(FPOP.models, file="FPOP.models.RData")
