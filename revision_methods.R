# Shared methods for the targeted extension; original diagnostics are unchanged.
library(stats)
library(utils)
source('diagnostics_core.R')
revision_angles <- list(
 canonical=as.matrix(expand.grid(a=c(0,pi/2),b=c(0,pi/2))),
 common_rotations=do.call(rbind,lapply(c(0,pi/8,pi/3,pi/4),function(a) as.matrix(expand.grid(a=c(0,pi/2),b=c(0,pi/2)))+a)),
 crossed_pi4=as.matrix(expand.grid(a=c(0,pi/4),b=c(0,pi/4))),
 crossed_band6=as.matrix(expand.grid(a=c(0,pi/7),b=c(0,pi/7))),
 crossed_irrational=as.matrix(expand.grid(a=c(0,pi/sqrt(2)),b=c(0,pi/sqrt(2)))))
revision_mmd <- function(x,y,z,signs) {
 matrices<-lapply(c(.3,.5,1),function(s) {
  kernel<-function(a,b) exp(-(1-cos(outer(a,b,'-')))/s^2)
  yz<-kernel(y,z)
  kernel(x,x)*(kernel(y,y)+kernel(z,z)-yz-t(yz))
 })
 matrices[[4]]<-Reduce('+',matrices)/3
 names(matrices)<-c('mmd03','mmd05','mmd10','mmd_mix')
 vapply(matrices,function(H) {
  observed<-sum(H); randomized<-colSums(signs*(H%*%signs))
  (1+sum(randomized>=observed-1e-11*max(1,abs(observed))))/(ncol(signs)+1)
 },numeric(1))
}
revision_field_single <- function(x,y,z,hx,ht) {
 gx<-2*pi*(0:15)/16;gy<-2*pi*(0:31)/32
 K<-wrapped_kernel(gx,x,hx); w<-K/rowSums(K); active<-which(1/rowSums(w^2)>=20)
 A<-wrapped_kernel(gy,z,ht); B<-wrapped_kernel(gy,y,ht)
 A<-sweep(A,2,colSums(A)*2*pi/length(gy),'/');B<-sweep(B,2,colSums(B)*2*pi/length(gy),'/')
 if(!length(active))return(matrix(0,1,length(x)))
 do.call(rbind,lapply(active,function(i)sweep(A-B,2,w[i,],'*')))
}
revision_evaluate <- function(x,y,z,signs,individual=FALSE) {
 p<-vapply(revision_angles,function(angles)sign_test(projection_features(x,y,z,angles),signs)$p,numeric(1))
 p<-c(p,fourier8=sign_test(fourier_features(x,y,z,8),signs)$p,
      fourier12=sign_test(fourier_features(x,y,z,12),signs)$p,
      field_max=sign_test(field_features(x,y,z)$C,signs)$p,revision_mmd(x,y,z,signs))
 if(individual) for(hx in c(.35,.7))for(ht in c(.15,.3)) {
  p[paste0('field_',hx,'_',ht)]<-sign_test(revision_field_single(x,y,z,hx,ht),signs)$p
 }
 p
}
revision_summary <- function(raw,groups) {
 z<-qnorm(.975)
 out<-lapply(split(raw,interaction(raw[groups],drop=TRUE)),function(d) {
  n<-nrow(d); rate<-mean(d$reject); den<-1+z^2/n
  mid<-(rate+z^2/(2*n))/den; half<-z*sqrt(rate*(1-rate)/n+z^2/(4*n^2))/den
  cbind(d[1,groups,drop=FALSE],data.frame(datasets=n,rejections=sum(d$reject),rate=rate,
   mcse=sqrt(rate*(1-rate)/n),lower=mid-half,upper=mid+half))
 })
 do.call(rbind,out)
}
revision_paired <- function(raw,groups) {
 out<-lapply(split(raw,interaction(raw[groups],drop=TRUE)),function(d) {
  base<-d[d$method=='canonical',c('rep','reject')];names(base)[2]<-'baseline'
  x<-merge(d,base,by='rep');x$delta<-as.numeric(x$reject)-as.numeric(x$baseline)
  do.call(rbind,lapply(split(x,x$method),function(v)cbind(v[1,c(groups,'method'),drop=FALSE],
    data.frame(datasets=nrow(v),difference=mean(v$delta),mcse=sd(v$delta)/sqrt(nrow(v))))))
 })
 do.call(rbind,out)
}
