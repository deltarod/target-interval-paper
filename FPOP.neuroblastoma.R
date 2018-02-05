source("packages.R")

data(neuroblastoma, package="neuroblastoma")
data(neuroblastomaProcessed, package="penaltyLearning")

dir.create("neuroblastoma")

f <- function(x)x[is.finite(x)]
min(f(neuroblastomaProcessed$target.mat[,1]))
max(f(neuroblastomaProcessed$target.mat[,2]))
penalty.grid <- exp(seq(-5, 7))

label.dt <- data.table(neuroblastoma$annotations)
profile.dt <- data.table(neuroblastoma$profiles)
setkey(profile.dt, profile.id, chromosome)

label.i <- label.dt[,which(profile.id==4 & chromosome==2)]

FPOP.neuroblastoma.list <- list()
for(label.i in 1:nrow(label.dt)){
  this.label <- label.dt[label.i]
  this.meta <- this.label[, list(profile.id, chromosome)]
  pid.chr <- this.meta[, paste0(profile.id, ".", chromosome)]
  cache.rds <- paste0("neuroblastoma/", pid.chr, ".rds")
  if(file.exists(cache.rds)){
    this.result <- readRDS(cache.rds)
  }else{
    cat(sprintf("%4d / %4d %s\n", label.i, nrow(label.dt), pid.chr))
    this.profile <- profile.dt[this.meta]
    data.vec <- this.profile$logratio
    label.df <- this.label
    getError <- function(penalty){
      stopifnot(is.numeric(penalty))
      stopifnot(length(penalty) == 1)
      result <- fpop::Fpop(data.vec, penalty)
      end.vec <- result$t.est
      change.vec <- end.vec[-length(end.vec)]
      rss.dt <- data.table(start=c(1, change.vec+1), end=end.vec)[, {
        pro <- data.vec[start:end]
        m <- mean(pro)
        res.vec <- pro-m
        list(rss=sum(res.vec * res.vec))
      }, by=list(start, end)]
      model.dt <- data.table(
        this.meta,
        penalty,
        total.loss=result$J.est,
        loss=sum(rss.dt$rss),
        n.changes=length(change.vec))
      change.dt <- if(length(change.vec)){
        data.table(
          model.dt,
          change=(this.profile[change.vec, position]+this.profile[change.vec+1, position])/2)
      }else{
        data.table(profile.id=character(), chromosome=character(), n.changes=numeric(), change=numeric())
      }
      L <- penaltyLearning::labelError(
        model.dt, this.label, change.dt,
        change.var="change",
        model.vars="n.changes",
        problem.vars=c("profile.id", "chromosome"))
      L$model.errors
    }
    error.list <- list()
    next.pen <- c(0, 1e100)
    iteration <- 0
    last.target.vec <- c(-Inf, Inf)
    target.result.list <- list()
    while(length(next.pen)){
      cat(
        "Next =", paste(next.pen, collapse=", "),
        "mc.cores=", getOption("mc.cores"),
        "\n")
      iteration <- iteration+1
      error.list[paste(next.pen)] <- mclapply.or.stop(next.pen, getError)
      error.dt <- do.call(rbind, error.list)[order(penalty)]
      print(error.dt[,.(penalty, loss, n.changes, fp, fn, errors)])
      unique.segs <- error.dt[, data.table(
        .SD[which.min(iteration)],
        penalties=.N
        ), by=list(n.changes)]
      path.dt <- data.table(penaltyLearning::modelSelection(
        unique.segs, "loss", "n.changes"))
      path.dt[, next.pen := max.lambda]
      path.dt[, already.computed := next.pen %in% names(error.list)]
      path.dt[, no.next := c(diff(n.changes) == -1, NA)]
      path.dt[, done := already.computed | no.next]
      path.dt[, is.min := errors==min(errors)]
      path.dt[, min.err.interval := cumsum(ifelse(
                                   c(is.min[1], diff(is.min))==1, 1, 0))]
      other.candidates <- path.dt[which(0<diff(fn) & diff(fp)<0)]
      interval.dt <- path.dt[is.min==TRUE, {
        i <- if(1 == .N || 0 == errors[1]){
          ## No middle candidate if there is only one model in the
          ## interval, or if there are no errors.
          NA
        }else{
          d <- data.table(
            i=1:(.N-1),
            ## do not attempt to explore other.candidates -- try
            ## different ones!
            is.other=next.pen[-.N] %in% other.candidates$next.pen,
            dist=diff(max.log.lambda)+diff(min.log.lambda),
            done=done[-.N])
          d[is.other==FALSE & done==FALSE, i[which.max(dist)]]
        }
        if(length(i)==0)i <- NA
        data.table(
          min.lambda=min.lambda[1],
          min.log.lambda=min.log.lambda[1],
          mid.lambda=max.lambda[i],
          max.lambda=max.lambda[.N],
          max.log.lambda=max.log.lambda[.N],
          log.size=max.log.lambda[.N]-min.log.lambda[1]
          )
      }, by=list(min.err.interval)]
      largest.interval <- interval.dt[which.max(log.size)]
      target.vec <- largest.interval[, c(min.log.lambda, max.log.lambda)]
      diff.target.vec <- target.vec-last.target.vec
      last.target.vec <- target.vec
      target.result.list[[paste(iteration)]] <- largest.interval[, data.table(
        iteration,
        min.log.lambda,
        max.log.lambda)]
      target.lambda <- largest.interval[, c(min.lambda, max.lambda)]
      error.candidates <- path.dt[next.pen %in% target.lambda]
      stopping.candidates <- rbind(error.candidates, other.candidates)[done==FALSE]
      cat(sprintf(
        "Target interval: %f %f change: %f %f\n",
        target.vec[1], target.vec[2],
        diff.target.vec[1], diff.target.vec[2]))
      next.pen <- if(nrow(stopping.candidates)){
        lambda.vec <- interval.dt[, c(min.lambda, mid.lambda, max.lambda)]
        interval.candidates <- path.dt[next.pen %in% lambda.vec][done==FALSE]
        unique(rbind(stopping.candidates, interval.candidates)$next.pen)
      }
    }#while(!is.null(pen))
    target.iterations <- do.call(rbind, target.result.list)
    grid.models <- do.call(rbind, lapply(penalty.grid, function(penalty){
      result <- fpop::Fpop(data.vec, penalty)
      end.vec <- result$t.est
      change.vec <- end.vec[-length(end.vec)]
      rss.dt <- data.table(start=c(1, change.vec+1), end=end.vec)[, {
        pro <- data.vec[start:end]
        m <- mean(pro)
        res.vec <- pro-m
        list(rss=sum(res.vec * res.vec))
      }, by=list(start, end)]
      model.dt <- data.table(
        this.meta,
        penalty,
        total.loss=result$J.est,
        loss=sum(rss.dt$rss),
        n.changes=length(change.vec))
      change.dt <- if(length(change.vec)){
        data.table(
          model.dt,
          change=(this.profile[change.vec, position]+this.profile[change.vec+1, position])/2)
      }else{
        data.table(profile.id=character(), chromosome=character(), n.changes=numeric(), change=numeric())
      }
      L <- penaltyLearning::labelError(
        model.dt, this.label, change.dt,
        change.var="change",
        model.vars="n.changes",
        problem.vars=c("profile.id", "chromosome"))
      L$model.errors
    }))
    grid.ms <- penaltyLearning::modelSelection(grid.models, "loss", "n.changes")
    grid.interval <- penaltyLearning::targetIntervals(grid.ms, problem.vars=c("profile.id", "chromosome"))
    PDPA.interval <- neuroblastomaProcessed$target.mat[pid.chr,]
    PDPA.errors <- neuroblastomaProcessed$errors[this.meta]
    this.result <- rbind(
      data.table(algo="grid.search", grid.interval, n.models=length(penalty.grid)),
      data.table(
        algo="PDPA", this.meta,
        min.log.lambda=PDPA.interval[1],
        max.log.lambda=PDPA.interval[2],
        errors=min(PDPA.errors$errors),
        n.models=20),
      data.table(
        algo="proposed", this.meta,
        largest.interval[, list(min.log.lambda, max.log.lambda)],
        errors=min(error.dt$errors), n.models=length(error.list)))
    saveRDS(this.result, cache.rds)
  }
  FPOP.neuroblastoma.list[[pid.chr]] <- this.result
}
FPOP.neuroblastoma <- do.call(rbind, FPOP.neuroblastoma.list)

save(FPOP.neuroblastoma, penalty.grid, file="FPOP.neuroblastoma.RData")
