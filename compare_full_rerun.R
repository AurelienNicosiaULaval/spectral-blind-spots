# Compare a complete recalculation (benchmark.R --output=DIR) with the recorded outcomes.
# Usage: Rscript compare_full_rerun.R DIR [OUTPUT_CSV]. Elapsed times are excluded.
args <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")][1])
root <- dirname(normalizePath(script))
new <- read.csv(file.path(args[1],"all_results.csv"))
old <- read.csv(file.path(root,"resultats_simulation","all_results.csv"))
key <- c("condition","replication","diagnostic")
stopifnot(nrow(new)==56000L,nrow(old)==56000L,!anyDuplicated(new[key]))
m <- merge(old,new,by=key,suffixes=c("_recorded","_rerun"))
stopifnot(nrow(m)==56000L)
by_diag <- do.call(rbind,lapply(split(m,m$diagnostic),function(g) data.frame(
  diagnostic=g$diagnostic[1],outcomes=nrow(g),
  seed_mismatches=sum(g$seed_recorded!=g$seed_rerun),
  decision_mismatches=sum(g$reject_recorded!=g$reject_rerun),
  p_value_mismatches=sum(g$p_value_recorded!=g$p_value_rerun),
  max_abs_p_difference=max(abs(g$p_value_recorded-g$p_value_rerun)),
  component_mismatches=sum(g$components_recorded!=g$components_rerun),
  max_relative_statistic_difference=max(abs(g$statistic_recorded-g$statistic_rerun)/pmax(1,abs(g$statistic_recorded))),
  field_position_mismatches=sum(g$field_min_positions_recorded!=g$field_min_positions_rerun |
                                g$field_max_positions_recorded!=g$field_max_positions_rerun))))
print(by_diag,row.names=FALSE)
out <- if(length(args)>1) args[2] else file.path(args[1],"full_rerun_comparison.csv")
write.csv(by_diag,out,row.names=FALSE)
cat(sprintf("Datasets: %d; outcomes: %d; decision mismatches: %d; p-value mismatches: %d.\n",
            length(unique(paste(m$condition,m$replication))),nrow(m),
            sum(by_diag$decision_mismatches),sum(by_diag$p_value_mismatches)))
