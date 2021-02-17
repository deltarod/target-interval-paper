source("packages.R")

u.vec <- c(
  "https://raw.githubusercontent.com/tdhock/feature-learning-benchmark/master/labeled_problems_folds.csv",
  "https://raw.githubusercontent.com/tdhock/feature-learning-benchmark/master/labeled_problems_features.csv",
  "https://raw.githubusercontent.com/tdhock/feature-learning-benchmark/master/labeled_problems_possible_errors.csv",
  "http://members.cbio.ensmp.fr/~thocking/data/target.intervals.models.csv")
data.list <- list()
for(u in u.vec){
  f <- basename(u)
  if(!file.exists(f)){
    download.file(u, f)
  }
  data.list[[f]] <- data.table::fread(f)
}
possible.dt <- nc::capture_first_df(
  data.list[["labeled_problems_possible_errors.csv"]],
  prob.dir=list(
    set.name="[^/]+",
    "/.*/",
    problem=".*"))
possible.problems.folds <- data.list[["labeled_problems_folds.csv"]][
  possible.dt, on=c("set.name", "problem")]

labeled_problems_penalties <- data.list[["target.intervals.models.csv"]][
  data.list[["labeled_problems_features.csv"]][["prob.dir"]], on="prob.dir"
][
  possible.problems.folds, on="prob.dir"
]

if(file.exists("labeled_problems_iterations.csv")){
  labeled_problems_iterations <- data.table::fread(
    "labeled_problems_iterations.csv")
}else{
  labeled_problems_iterations <- labeled_problems_penalties[, {
    one.prob <- data.table(iteration, peaks, total.cost, errors)
    one.ms <- data.table(max.it=unique(iteration))[, {
      some.models <- one.prob[iteration <= max.it]
      penaltyLearning::modelSelection(some.models, "total.cost", "peaks")
    }, by=max.it]
    it.dt <- penaltyLearning::targetIntervals(one.ms, "max.it")[order(max.it)]
    for(prefix in c("min", "max")){
      limit <- it.dt[[paste0(prefix, ".log.lambda")]]
      d <- ifelse(limit[-1]==limit[-length(limit)], 0, abs(diff(limit)))
      set(
        it.dt,
        j=paste0(prefix, ".change"),
        value=c(NA, d))
    }
    it.dt[, total.change := min.change + max.change]
    it.dt[, .(max.it, errors, total.change)]
  }, by=prob.dir]
  data.table::fwrite(
    labeled_problems_iterations,
    "labeled_problems_iterations.csv")
}
labeled_problems_iterations[, cummin := c(NA, cummin(total.change[-1])), by=prob.dir]

ggplot()+
  geom_line(aes(
    max.it, total.change, group=prob.dir),
    data=labeled_problems_iterations)+
  scale_y_log10()

ggplot()+
  geom_line(aes(
    max.it, errors, group=prob.dir),
    data=labeled_problems_iterations)+
  scale_y_log10()

## strategy: keep doing iterations until cummin drops below thresh.
thresh.vec <- c(0, 0.001, 0.01, 0.1, 0.3, 0.5, 1, 1.5, 2, Inf)
it.per.prob <- data.table(thresh=thresh.vec)[, {
  labeled_problems_iterations[, {
    below.i <- which(cummin < thresh)
    .(max.it=if(length(below.i))below.i[1] else .N)
  }, by=prob.dir]
}, by=thresh]
med.dt <- it.per.prob[, .(
  median=median(max.it),
  mean=mean(max.it)
), by=thresh]
med.color <- "red"
gg <- ggplot()+
  geom_bar(aes(
    max.it),
    data=it.per.prob)+
  geom_vline(aes(
    xintercept=mean),
    color=med.color,
    data=med.dt)+
  geom_text(aes(
    mean, 0, label=sprintf(" mean=%.1f", mean)),
    color=med.color,
    vjust=0,
    hjust=0,
    data=med.dt)+
  facet_grid(thresh ~ ., labeller=label_both, scales="free")
png("figure-penalties-per-problem-maxit-thresh-hist.png", width=13, height=8, units="in", res=100)
print(gg)
dev.off()

## Do error rates always decrease after the first iteration? no.
err.ranges <- labeled_problems_iterations[, .(
  first.errors=errors[1],
  min.errors=min(errors)
), by=prob.dir]
err.ranges[first.errors==min.errors]

possible.folds <- possible.problems.folds[, .(
  possible.fp=sum(possible.fp),
  possible.tp=sum(possible.tp)
), by=c("set.name", "fold")]

penalties.dt.list <- list()
for(test.fold.i in 1:nrow(possible.folds)){
  cat(sprintf("%4d / %4d folds\n", test.fold.i, nrow(possible.folds)))
  test.fold.row <- possible.folds[test.fold.i]
  set.penalties <- labeled_problems_penalties[
    set.name == test.fold.row$set.name]
  train.penalties <- set.penalties[
    fold != test.fold.row$fold]
  train.ord <- train.penalties[order(
    iteration, # go in order of target interval search algo.
    -penalty, # break ties with 0/Inf-> Inf first.
    bedGraph.lines # smaller contigs first.
  )]
  train.ord[, cum.seconds := cumsum(seconds)]
  prob.seconds <- train.ord[, .(
    sum.seconds=sum(seconds)
  ), by=.(prob.dir, bedGraph.lines)][order(bedGraph.lines)]
  all.n <- 1:nrow(prob.seconds)
  target.n <- c(
    1:9,
    10^seq(1, log10(nrow(prob.seconds)), by=0.5),
    nrow(prob.seconds))
  some.n <- target.n[target.n %in% all.n]
  for(n.full.problems in some.n){
    full.probs <- prob.seconds[1:n.full.problems]
    n.probs.seconds <- sum(full.probs$sum.seconds)
    penalties.list <- list(
      full=train.penalties[full.probs$prob.dir, on="prob.dir"],
      partial=train.ord[cum.seconds <= n.probs.seconds])
    for(order.type in names(penalties.list)){
      order.dt <- penalties.list[[order.type]]
      prob.counts <- order.dt[, .(
        penalties=.N,
        min.errors=min(errors),
        max.errors=max(errors)
      ), by=prob.dir]
      n.train <- prob.counts[, sum(min.errors < max.errors)]
      range.dt <- prob.counts[, .(
        min=min(penalties),
        max=max(penalties)
      )]
      n.problems <- nrow(prob.counts)
      n.penalties <- nrow(order.dt)
      ## TODO train on order.dt, compute test accuracy/AUC.
      penalties.dt.list[[paste(test.fold.i, n.full.problems, order.type)]] <- data.table(
        test.fold.row,
        n.full.problems,
        order.type,
        n.penalties,
        n.train,
        range.dt)
    }
  }
}
(penalties.dt <- do.call(rbind, penalties.dt.list))

library(ggplot2)
gg <- ggplot()+
  geom_line(aes(
    n.full.problems, n.train,
    group=paste(order.type, fold),
    color=order.type),
    data=penalties.dt)+
  facet_wrap("set.name", scales="free")+
  scale_x_log10()+
  scale_y_log10()
png("figure-penalties-per-problem-train.png", width=13, height=8, units="in", res=100)
print(gg)
dev.off()

gg <- ggplot()+
  geom_line(aes(
    n.full.problems, n.penalties,
    group=paste(order.type, fold),
    color=order.type),
    data=penalties.dt)+
  facet_wrap("set.name", scales="free")+
  scale_x_log10()+
  scale_y_log10()
png("figure-penalties-per-problem.png", width=13, height=8, units="in", res=100)
print(gg)
dev.off()
