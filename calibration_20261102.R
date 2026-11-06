library(stats)
library(utils)
library(parallel)
source('diagnostics_core.R')
RNGkind('Mersenne-Twister','Inversion','Rejection')
dir.create('renforcement',showWarnings=FALSE)
draw_vm_revision<-function(x) {
 out<-numeric(length(x));pending<-seq_along(x)
 while(length(pending)) {
  y<-runif(length(pending),0,2*pi);rad<-sqrt(4+sin(x[pending])^2)
  keep<-runif(length(pending))<exp(2*cos(y)+sin(x[pending])*sin(y)-rad)
  out[pending[keep]]<-y[keep];pending<-pending[!keep]
 }
 out
}
rows<-mclapply(1:2000,function(rep) {
 seed<-2026211000L+rep;set.seed(seed);x<-runif(240,0,2*pi)
 y<-draw_vm_revision(x);zz<-draw_vm_revision(x)
 signs<-matrix(sample(c(-1,1),240*199,TRUE),240,199)
 p<-sign_test(fourier_features(x,y,zz,12),signs)$p
 data.frame(rep=rep,seed=seed,p=p,reject=p<=.05)
},mc.cores=if(.Platform$OS.type=='windows')1L else 2L,mc.set.seed=FALSE)
raw<-do.call(rbind,rows);stopifnot(nrow(raw)==2000,all(is.finite(raw$p)))
write.csv(raw,'renforcement/calibration_2000.csv',row.names=FALSE)
prior<-read.csv('complement/resultats.csv');prior<-prior[prior$case=='fidelity'&prior$method=='fourier12',]
extra<-read.csv('complement/calibration_supplementaire.csv')
summary<-data.frame(run=c('initial_post_review','independent_500','new_2000','pooled'),
 n=c(nrow(prior),nrow(extra),nrow(raw),nrow(prior)+nrow(extra)+nrow(raw)),
 rejections=c(sum(prior$reject),sum(extra$reject),sum(raw$reject),sum(prior$reject)+sum(extra$reject)+sum(raw$reject)))
z<-qnorm(.975);p<-summary$rejections/summary$n;den<-1+z^2/summary$n
mid<-(p+z^2/(2*summary$n))/den;half<-z*sqrt(p*(1-p)/summary$n+z^2/(4*summary$n^2))/den
summary$rate<-p;summary$lower<-mid-half;summary$upper<-mid+half
summary$upper_tail<-mapply(function(k,n)pbinom(k-1,n,.05,lower.tail=FALSE),summary$rejections,summary$n)
write.csv(summary,'renforcement/calibration_synthese.csv',row.names=FALSE)
print(summary)
