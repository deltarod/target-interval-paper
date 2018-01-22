library(data.table)
library(PeakError)
counts.RData.vec <- Sys.glob("data/H*/*/counts.RData")
PDPA.models.list <- list()
i.vec <- seq_along(counts.RData.vec)
i.vec <- 1:2
PDPA.modelSelection.dt.list <- list()
PDPA.targets.dt.list <- list()
for(counts.RData.i in i.vec){
  counts.RData <- counts.RData.vec[[counts.RData.i]]
  chunk.dir <- dirname(counts.RData)
  load(counts.RData)
  counts.dt <- data.table(counts)
  setkey(counts.dt, sample.id)
  chunk.name <- sub("data/", "", chunk.dir)
  set.name <- dirname(chunk.name)
  chunk.id <- basename(chunk.name)
  PDPA.model.RData <- file.path(chunk.dir, "PDPA.model.RData")
  if(file.exists(PDPA.model.RData)){
    load(PDPA.model.RData)
    if("count.vec" %in% names(PDPA.model[[1]])){
      unlink(PDPA.model.RData)
    }
  }
  if(!file.exists(PDPA.model.RData)){
    u <- sub(
      "data/",
      "http://members.cbio.mines-paristech.fr/~thocking/chip-seq-chunk-db/",
      PDPA.model.RData)
    download.file(u, PDPA.model.RData)
  }
  load(PDPA.model.RData)
  regions.RData <- file.path(chunk.dir, "regions.RData")
  load(regions.RData)
  regions.dt <- data.table(regions)
  setkey(regions.dt, sample.id)
  PDPA.errors.list <- list()
  for(sample.id in names(PDPA.model)){
    sample.model <- PDPA.model[[sample.id]]
    sample.counts <- counts.dt[sample.id]
    sample.regions <- regions.dt[sample.id]
    for(n.segments in seq(1, 19, by=2)){
      change.vec <- sample.model$ends.mat[n.segments,2:n.segments]
      end.vec <- c(change.vec, nrow(sample.counts))
      start.vec <- c(1, change.vec+1)
      seg.dt <- sample.counts[, data.table(
        chromStart=chromStart[start.vec],
        chromEnd=chromEnd[end.vec],
        status=rep(c("background", "peak"), l=n.segments))]
      peak.dt <- seg.dt[status=="peak"]
      error.df <- PeakErrorChrom(peak.dt, sample.regions)
      PDPA.errors.list[[paste(sample.id, n.segments)]] <- with(error.df, data.table(sample.id, 
                                                                                    peaks=(n.segments-1)/2,
                                                                                    loss=sample.model$loss.vec[n.segments],
        errors=sum(fp+fn)))
    }
  }
  PDPA.errors <- do.call(rbind, PDPA.errors.list)
  PDPA.modelSelection <- PDPA.errors[, penaltyLearning::modelSelection(
    .SD, "loss", "peaks"), by=sample.id]
  PDPA.modelSelection.dt.list[[paste(chunk.dir)]] <- data.table(
    chunk.name, set.name, chunk.id, PDPA.modelSelection)
  PDPA.targets <- penaltyLearning::targetIntervals(
    PDPA.modelSelection, "sample.id")
  PDPA.targets.dt.list[[paste(chunk.dir)]] <- data.table(
    chunk.name, set.name, chunk.id, PDPA.targets)
}
PDPA.targets.dt <- do.call(rbind, PDPA.targets.dt.list)
PDPA.modelSelection.dt <- do.call(rbind, PDPA.modelSelection.dt.list)

save(PDPA.targets.dt, PDPA.modelSelection.dt, file="PDPA.models.RData")
