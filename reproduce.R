# Rebuild summaries and figures from recorded outcomes.
library(utils)

arguments <- commandArgs(trailingOnly = TRUE)
allowed <- c("--verify", "--help")
if (any(!arguments %in% allowed)) {
  stop("Allowed options: --verify, --help", call. = FALSE)
}
if ("--help" %in% arguments) {
  cat("Usage: Rscript --vanilla reproduce.R [--verify]\n",
      "Rebuilds checks, examples, tables and figures.\n",
      "--verify also recalculates 44 datasets and compares 440 outcomes.\n",
      "See README.md for the full simulation.\n", sep = "")
  quit(save = "no", status = 0L)
}

script_argument <- commandArgs()[startsWith(commandArgs(), "--file=")]
stopifnot(length(script_argument) == 1L)
script_path <- sub("^--file=", "", script_argument)
# R encodes spaces in --file arguments as ~+~ on some platforms.
script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
root <- dirname(normalizePath(script_path))
setwd(root)
rscript <- file.path(R.home("bin"), "Rscript")
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 before rebuilding the figures.", call. = FALSE)
}

run_step <- function(script, options = character()) {
  cat("\nRunning ", script, "\n", sep = "")
  status <- system2(rscript, c("--vanilla", shQuote(script), options))
  if (status != 0L) {
    stop(paste(script, "failed with exit code", status), call. = FALSE)
  }
}

steps <- c("check_geometry.R", "check_nullspace.R", "check_diagnostics.R",
           "make_example.R", "toy.R", "report_results.R")
for (script in steps) run_step(script)
if ("--verify" %in% arguments) run_step("benchmark.R", "--verify")

cat("\nTables and figures rebuilt.\n")
