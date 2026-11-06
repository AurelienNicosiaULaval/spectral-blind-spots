library(stats)
library(utils)
library(parallel)
source('revision_methods.R')
RNGkind('Mersenne-Twister','Inversion','Rejection')
outdir<-'renforcement';dir.create(outdir,showWarnings=FALSE)
cases<-expand.grid(model=c('wc05','wc08','wc09','axial08'),fidelity=c(TRUE,FALSE),stringsAsFactors=FALSE)
cases$rho<-c(.5,.8,.9,.8)[match(cases$model,c('wc05','wc08','wc09','axial08'))]
cases$condition<-seq_len(nrow(cases));cases$delta<-.5
write.csv(cases,file.path(outdir,'rotation_conditions.csv'),row.names=FALSE)
# Inverse CDF sampling; the axial version adds an independent half-turn.
draw_join <- function(n,rho,shift,axial=FALSE) {
 u<-runif(n,-pi,pi)
 t<-2*atan((1-rho)/(1+rho)*tan(u/2))
 if(axial)t<-t/2+sample(c(0,pi),n,TRUE)
 wrap(t+shift)
}
run_condition<-function(i) {
 ca<-cases[i,];file<-file.path(outdir,paste0('rotation_checkpoint_',i,'.rds'))
 if(file.exists(file))return(readRDS(file))
 rows<-lapply(seq_len(500),function(rep) {
  seed<-2026110200L+10000L*i+rep;set.seed(seed)
  x<-runif(240,0,2*pi);y<-wrap(x+draw_join(240,ca$rho,.5,ca$model=='axial08'))
  shift<-if(ca$fidelity).5 else -.5
  zz<-wrap(x+draw_join(240,ca$rho,shift,ca$model=='axial08'))
  signs<-matrix(sample(c(-1,1),240*199,TRUE),240,199)
  p<-revision_evaluate(x,y,zz,signs)
  stopifnot(length(p)==12,all(is.finite(p)),all(p>=.005&p<=1))
  data.frame(condition=i,model=ca$model,fidelity=ca$fidelity,rho=ca$rho,delta=.5,
             n=240,rep=rep,seed=seed,method=names(p),p=unname(p),reject=unname(p)<=.05)
 })
 result<-do.call(rbind,rows);saveRDS(result,file);result
}
raw<-do.call(rbind,mclapply(seq_len(nrow(cases)),run_condition,mc.cores=if(.Platform$OS.type=='windows')1L else 2L,mc.set.seed=FALSE))
stopifnot(nrow(raw)==48000,!anyNA(raw))
write.csv(raw,file.path(outdir,'rotation_resultats.csv'),row.names=FALSE)
write.csv(revision_summary(raw,c('condition','model','fidelity','method')),file.path(outdir,'rotation_synthese.csv'),row.names=FALSE)
write.csv(revision_paired(raw,c('condition','model','fidelity')),file.path(outdir,'rotation_ecarts_apparies.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(outdir,'rotation_sessionInfo.txt'))
cat('Rotation study: 4000 datasets, 48000 outcomes completed.\n')
