library(stats)
library(utils)
library(parallel)
source('revision_methods.R')
RNGkind('Mersenne-Twister','Inversion','Rejection')
dir.create('renforcement',showWarnings=FALSE)
recorded<-read.csv('resultats_simulation/all_results.csv')
cases<-unique(recorded[recorded$scenario%in%c('phase_4','phase_6'),c('scenario','condition','n')])
rows<-mclapply(seq_len(nrow(cases)),function(i) {
 ca<-cases[i,]; checkpoint<-paste0('renforcement/sensibilite_checkpoint_',ca$condition,'.rds')
 if(file.exists(checkpoint))return(readRDS(checkpoint))
 res<-lapply(1:200,function(rep) {
  old<-recorded[recorded$condition==ca$condition & recorded$replication==rep,]
  seed<-unique(old$seed);stopifnot(length(seed)==1);set.seed(seed)
  x<-draw_x(ca$n,0);definition<-list(kind='phase',order=as.integer(sub('phase_','',ca$scenario)),amplitude=.8)
  y<-draw_response(x,definition);zz<-runif(ca$n,0,2*pi)
  signs<-matrix(sample(c(-1,1),ca$n*199,TRUE),ca$n,199)
  p<-revision_evaluate(x,y,zz,signs,individual=TRUE)
  map<-c(canonical='canonical',common_rotations='common_rotations',crossed_pi4='crossed_pi4',crossed_band6='crossed_band6',crossed_irrational='crossed_irrational',fourier8='fourier_8',field_max='periodic_field',mmd05='paired_mmd')
  expected<-old$p_value[match(unname(map),old$diagnostic)]
  stopifnot(max(abs(p[names(map)]-expected))<1e-12)
  data.frame(scenario=ca$scenario,condition=ca$condition,n=ca$n,rep=rep,seed=seed,
   method=names(p),p=unname(p),reject=unname(p)<=.05,replay_error=max(abs(p[names(map)]-expected)))
 })
 result<-do.call(rbind,res);saveRDS(result,checkpoint);result
},mc.cores=if(.Platform$OS.type=='windows')1L else 2L,mc.set.seed=FALSE)
raw<-do.call(rbind,rows);stopifnot(nrow(raw)==12800,all(is.finite(raw$p)))
write.csv(raw,'renforcement/sensibilite_resultats.csv',row.names=FALSE)
write.csv(revision_summary(raw,c('scenario','n','method')),'renforcement/sensibilite_synthese.csv',row.names=FALSE)
write.csv(revision_paired(raw,c('scenario','n')),'renforcement/sensibilite_ecarts_apparies.csv',row.names=FALSE)
cat('Sensitivity: 800 original datasets replayed, 6400 original p-values matched, 12800 extended outcomes.\n')
