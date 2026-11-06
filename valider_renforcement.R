library(stats)
library(utils)
source('revision_methods.R')
rot<-read.csv('renforcement/rotation_resultats.csv')
sen<-read.csv('renforcement/sensibilite_resultats.csv')
cal<-read.csv('renforcement/calibration_2000.csv')
stopifnot(nrow(rot)==48000,nrow(sen)==12800,nrow(cal)==2000,
          length(unique(paste(rot$condition,rot$rep)))==4000,
          length(unique(paste(sen$condition,sen$rep)))==800,
          length(intersect(unique(rot$seed),cal$seed))==0,
          all(is.finite(rot$p)),all(is.finite(sen$p)),all(is.finite(cal$p)),
          max(sen$replay_error)<1e-12)
for(name in c('rotation','sensibilite')) {
 raw<-if(name=='rotation')rot else sen
 groups<-if(name=='rotation')c('condition','model','fidelity','method') else c('scenario','n','method')
 expected<-revision_summary(raw,groups)
 stored<-read.csv(paste0('renforcement/',name,'_synthese.csv'))
 expected<-expected[order(do.call(paste,expected[groups])),];stored<-stored[order(do.call(paste,stored[groups])),]
 stopifnot(isTRUE(all.equal(expected,stored,check.attributes=FALSE,tolerance=1e-12)))
}
original_seeds<-unique(read.csv('resultats_simulation/all_results.csv')$seed)
prior_calibration_seeds<-c(read.csv('complement/calibration_supplementaire.csv')$seed,
 subset(read.csv('complement/resultats.csv'),case=='fidelity' & method=='fourier12')$seed)
stopifnot(length(intersect(cal$seed,prior_calibration_seeds))==0,
          length(intersect(unique(rot$seed),original_seeds))==0)
energy<-read.csv('renforcement/rotation_energie.csv')
stopifnot(max(abs(energy$series-energy$grid256))<1e-8,
          max(abs(energy$series-energy$grid512))<1e-8,
          max(abs(energy$integral-1))<1e-10,max(energy$margin_error)<1e-10)
# Verify kernel-mixture exchange equivariance without relying on the simulated rates.
set.seed(2026110299);n<-20L;x<-runif(n,0,2*pi);y<-runif(n,0,2*pi);z<-runif(n,0,2*pi)
epsilon<-sample(c(-1,1),n,TRUE)
make_H<-function(y,z) {
 Reduce('+',lapply(c(.3,.5,1),function(s) {
  k<-function(a,b)exp(-(1-cos(outer(a,b,'-')))/s^2)
  k(x,x)*(k(y,y)+k(z,z)-k(y,z)-k(z,y))
 }))/3
}
H<-make_H(y,z);HS<-make_H(ifelse(epsilon==1,y,z),ifelse(epsilon==1,z,y))
stopifnot(max(abs(HS-H*outer(epsilon,epsilon)))<1e-12,
          min(eigen(H,symmetric=TRUE,only.values=TRUE)$values)>-1e-10)
# Conditional sampling checks, separately from the simulation outcome tables.
expr<-parse('etude_rotation_20261102.R')
for(e in expr)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('draw_join')))eval(e)
for(axial in c(FALSE,TRUE)) {
 z<-draw_join(50000,.8,.5,axial);frequency<-if(axial)2 else 1
 target<-.8*exp(1i*frequency*.5);values<-exp(1i*frequency*z)
 stopifnot(abs(mean(Re(values))-Re(target))<5*sd(Re(values))/sqrt(length(z)),
           abs(mean(Im(values))-Im(target))<5*sd(Im(values))/sqrt(length(z)))
}
scores<-read.csv('application/scores_paires.csv')
daily<-read.csv('application/scores_journaliers.csv')
stopifnot(nrow(scores)==276,nrow(daily)==23,
          max(abs(scores$difference-scores$log_persistence+scores$log_original))<1e-12,
          abs(mean(daily$difference)-mean(scores$difference))<1e-12)
v<-read.csv('application/vues_projetees.csv')
stopifnot(nrow(v)==8*24*24,max(abs(aggregate(v$difference,list(v$model,v$view),sum)$x))<1e-12)
best<-read.csv('application/bandes_cv.csv')
stopifnot(nrow(best)==25,sum(best$selected)==1,
          all(is.finite(best$loglik)),max(best$loglik)==best$loglik[best$selected])
original<-read.csv('application/tirages.csv')
new<-read.csv('application/tirages_persistance.csv')
stopifnot(identical(dim(original),c(276L,50L)),identical(dim(new),dim(original)))
writeLines(c('48000 rotation outcomes, 12800 sensitivity outcomes and 2000 new calibration outcomes verified.',
 'Series/quadrature, margins, exchange equivariance, mixture positivity and conditional draws checked.',
 'Wind scores, 25 CV scores, probability-bin totals and training/evaluation sizes checked.'),
 'renforcement/validation.txt')
cat(readLines('renforcement/validation.txt'),sep='\n')
