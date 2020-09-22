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

possible.folds <- possible.problems.folds[, .(
  possible.fp=sum(possible.fp),
  possible.tp=sum(possible.tp)
), by=c("set.name", "fold")]

penalties.dt.list <- list()
for(test.fold.i in 1:nrow(possible.folds)){
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
  for(n.problems in 1:nrow(prob.seconds)){
    full.probs <- prob.seconds[1:n.problems]
    n.probs.seconds <- sum(full.probs$sum.seconds)
    penalties.list <- list(
      full=train.penalties[full.probs$prob.dir, on="prob.dir"],
      partial=train.ord[cum.seconds <= n.probs.seconds])
    for(order.type in names(penalties.list)){
      order.dt <- penalties.list[[order.type]]
      ## TODO train on order.dt, compute test accuracy/AUC.
      penalties.dt.list[[paste(test.fold.i, n.problems, order.type)]] <- data.table(
        test.fold.row,
        n.problems,
        order.type,
        penalties=nrow(order.dt))
    }
  }
}
penalties.dt <- do.call(rbind, penalties.dt.list)

library(ggplot2)
gg <- ggplot()+
  geom_line(aes(
    n.problems, penalties,
    group=paste(order.type, fold),
    color=order.type),
    data=penalties.dt)+
  facet_wrap("set.name", scales="free")+
  scale_x_log10()+
  scale_y_log10()
png("figure-penalties-per-problem.png", width=13, height=8, units="in", res=100)
print(gg)
dev.off()
