library(stats)
library(utils)
library(parallel)
library(tools)
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")][1])
root <- dirname(normalizePath(script))
source(file.path(root,"diagnostics_core.R"))
arguments <- commandArgs(trailingOnly=TRUE)
output_argument <- arguments[startsWith(arguments,"--output=")]
stopifnot(length(output_argument)<=1L,
          all(arguments %in% c("--smoke","--verify") | startsWith(arguments,"--output=")))
out <- if(length(output_argument)) sub("^--output=","",output_argument) else
  file.path(root,"resultats_simulation")
dir.create(out,showWarnings=FALSE)

definitions <- data.frame(
  scenario=c("null_uniform","null_conditional","odd_1","phase_2_moderate",
             "phase_2_strong","phase_4","phase_6","phase_7_outside",
             "product_2","local_rare","local_common"),
  kind=c("null_uniform","null_conditional",rep("phase",6),"product",
         "local_rare","local_common"),
  order=c(0,1,1,2,2,4,6,7,2,2,2),
  amplitude=c(0,0.8,0.4,0.4,0.8,0.8,0.8,0.8,0.8,0.8,0.8),
  rho=c(0,0.75,rep(0,7),0.75,0.75),
  is_null=c(TRUE,TRUE,rep(FALSE,9))
)
conditions <- merge(definitions,data.frame(n=c(240L,720L)),sort=FALSE)
conditions <- conditions[order(match(conditions$scenario,definitions$scenario),conditions$n),]
rownames(conditions) <- NULL
conditions$condition <- seq_len(nrow(conditions))
conditions$replications <- ifelse(conditions$is_null,500L,200L)
stopifnot(sum(conditions$replications)==5600L)

crossed <- function(delta) as.matrix(expand.grid(a=c(0,delta),b=c(0,delta)))
angle_sets <- list(canonical=crossed(pi/2),
  common_rotations=do.call(rbind,lapply(c(0,pi/8,pi/3,pi/4),function(a) crossed(pi/2)+a)),
  crossed_pi4=crossed(pi/4),crossed_band6=crossed(pi/7),
  crossed_irrational=crossed(pi/sqrt(2)))

evaluate <- function(definition,seed) {
  set.seed(seed)
  n <- definition$n
  x <- draw_x(n,definition$rho)
  y <- draw_response(x,definition)
  z <- if(definition$is_null) draw_response(x,definition) else runif(n,0,2*pi)
  signs <- matrix(sample(c(-1,1),n*199,TRUE),n,199)
  rows <- list()
  record <- function(name,operation,view_count=NA_integer_) {
    start <- proc.time()[[3]]
    result <- operation()
    elapsed <- proc.time()[[3]]-start
    rows[[length(rows)+1]] <<- data.frame(diagnostic=name,p_value=result$p,
      statistic=result$statistic,components=result$components,views=view_count,
      elapsed_seconds=elapsed,reject=result$p<=0.05)
  }
  for(name in names(angle_sets)) {
    angles <- angle_sets[[name]]
    record(name,function() sign_test(projection_features(x,y,z,angles),signs),nrow(angles))
  }
  record("fourier_8",function() sign_test(fourier_features(x,y,z),signs))
  record("marginal_8",function() {
    values <- exp(1i*outer(1:8,z))-exp(1i*outer(1:8,y))
    sign_test(rbind(Re(values),Im(values)),signs)
  })
  field <- NULL
  record("periodic_field",function() {
    field <<- field_features(x,y,z)
    sign_test(field$C,signs)
  })
  stopifnot(field$identity_error<1e-12)
  record("paired_mmd",function() mmd_test(x,y,z,signs))
  record("coverage_known_model",function() coverage_test(x,y,definition))
  result <- do.call(rbind,rows)
  result$seed <- seed
  result$field_min_positions <- field$min_support
  result$field_max_positions <- field$max_support
  result$field_identity_error <- field$identity_error
  result
}

if("--smoke" %in% arguments) {
  rows <- lapply(c(1L,3L,9L),function(index) {
    definition <- conditions[index,]
    result <- evaluate(definition,48100L+index)
    result$scenario <- definition$scenario; result$n <- definition$n
    result
  })
  write.csv(do.call(rbind,rows),file.path(root,"smoke_development.csv"),row.names=FALSE)
  cat("Three development runs complete, excluded from the benchmark.\n")
  quit(save="no",status=0)
}

# Recalculate the first and last replication of every condition. Timing is excluded.
if("--verify" %in% arguments) {
  recorded <- read.csv(file.path(out,"all_results.csv"))
  errors <- list()
  for(index in seq_len(nrow(conditions))) {
    definition <- conditions[index,]
    for(replication in c(1L,definition$replications)) {
      seed <- 2026101000L+10000L*index+replication
      recalculated <- evaluate(definition,seed)
      expected <- recorded[recorded$condition==index & recorded$replication==replication,]
      expected <- expected[match(recalculated$diagnostic,expected$diagnostic),]
      stopifnot(nrow(expected)==10L,all(recalculated$seed==expected$seed),
                all(recalculated$reject==expected$reject),
                all(recalculated$components==expected$components))
      p_error <- max(abs(recalculated$p_value-expected$p_value))
      statistic_error <- max(abs(recalculated$statistic-expected$statistic)/
                               pmax(1,abs(expected$statistic)))
      stopifnot(p_error<1e-12,statistic_error<1e-12)
      errors[[length(errors)+1L]] <- data.frame(index,replication,seed,p_error,statistic_error)
    }
  }
  write.csv(do.call(rbind,errors),file.path(out,"replay_checks.csv"),row.names=FALSE)
  cat("44 datasets and 440 diagnostic outcomes recalculated and matched; timings excluded.\n")
  quit(save="no",status=0)
}

fingerprint <- unname(md5sum(file.path(root,"diagnostics_core.R")))
run_condition <- function(index) {
  definition <- conditions[index,]
  path <- file.path(out,sprintf("condition_%02d.rds",index))
  saved <- if(file.exists(path)) readRDS(path) else
    list(definition=definition,core_md5=fingerprint,results=NULL)
  stopifnot(identical(saved$definition,definition),identical(saved$core_md5,fingerprint))
  complete <- if(is.null(saved$results)) 0L else max(saved$results$replication)
  if(complete<definition$replications) for(replication in (complete+1L):definition$replications) {
    seed <- 2026101000L+10000L*index+replication
    result <- evaluate(definition,seed)
    result$condition <- index; result$scenario <- definition$scenario
    result$n <- definition$n; result$replication <- replication
    saved$results <- rbind(saved$results,result)
    if(replication %% 20L==0L || replication==definition$replications) {
      temporary <- paste0(path,".tmp")
      saveRDS(saved,temporary)
      stopifnot(file.rename(temporary,path))
    }
    if(replication %% 100L==0L || replication==definition$replications) {
      cat(sprintf("Condition %02d: %s, n=%d, %d/%d\n",index,definition$scenario,
                  definition$n,replication,definition$replications))
      flush.console()
    }
  }
  saved$results
}

write.csv(conditions,file.path(out,"conditions.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
if(!file.exists(file.path(out,"run_info.rds"))) saveRDS(list(started_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
             core_md5=fingerprint,conditions=conditions,signs=199L,alpha=0.05),
        file.path(out,"run_info.rds"))
workers <- if(.Platform$OS.type=="windows") 1L else 2L
results <- mclapply(seq_len(nrow(conditions)),run_condition,mc.cores=workers,mc.set.seed=FALSE)
failed <- vapply(results,inherits,logical(1),"try-error")
if(any(failed)) stop(paste("Conditions failed:",paste(which(failed),collapse=",")))
results <- do.call(rbind,results)
stopifnot(nrow(results)==56000L)
write.csv(results,file.path(out,"all_results.csv"),row.names=FALSE)
saveRDS(list(completed_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
             datasets=5600L,rows=nrow(results),core_md5=fingerprint),
        file.path(out,"completion.rds"))
cat("Benchmark completed: 5600 datasets, ten diagnostics.\n")
