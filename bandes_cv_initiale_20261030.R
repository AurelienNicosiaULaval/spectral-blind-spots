library(stats)
library(utils)
source('diagnostics_core.R')
p<-read.csv('application/paires.csv');days<-sort(unique(p$date));train<-p[p$date %in% days[1:50],];test<-p[p$date %in% days[51:73],];z<-as.matrix(read.csv('application/tirages.csv'))
stopifnot(nrow(train)==600,nrow(test)==276,ncol(z)==50)
cv<-expand.grid(hx=c(.35,.70,1.20),ht=c(.15,.30,.60))
cv$loglik<-vapply(seq_len(nrow(cv)),function(j){
 X<-wrapped_kernel(train$x,train$x,cv$hx[j]);Y<-wrapped_kernel(train$y,train$y,cv$ht[j]);X[outer(train$date,train$date,'==')]<-0
 density<-rowSums(X*Y)/rowSums(X);if(any(!is.finite(density)|density<=0)) stop('Invalid CV density')
 sum(log(density))
},0.0)
cv$selected<-seq_len(nrow(cv))==which.max(cv$loglik);write.csv(cv,'application/bandes_cv.csv',row.names=FALSE)
best<-cv[cv$selected,];gx<-2*pi*(0:31)/32;gy<-2*pi*(0:63)/64
field<-function(hx,ht){
 K<-wrapped_kernel(gx,test$x,hx);w<-K/rowSums(K);neff<-1/rowSums(w^2)
 ky<-wrapped_kernel(gy,test$y,ht);ky<-sweep(ky,2,colSums(ky)*2*pi/length(gy),'/')
 f<-matrix(0,length(gx),length(gy))
 for(j in seq_len(ncol(z))){kz<-wrapped_kernel(gy,z[,j],ht);kz<-sweep(kz,2,colSums(kz)*2*pi/length(gy),'/');f<-f+w%*%t(kz-ky)/ncol(z)}
 list(field=f,neff=neff)
}
a<-field(.7,.3);b<-field(best$hx,best$ht);mask<-a$neff>=20&b$neff>=20
if(!any(mask))stop('No shared support')
summary<-data.frame(hx=best$hx,ht=best$ht,loglik=best$loglik,active=sum(mask),minimum_support=min(b$neff),rms_diff=sqrt(mean((a$field[mask,]-b$field[mask,])^2)),max_abs_diff=max(abs(a$field[mask,]-b$field[mask,])))
write.csv(summary,'application/bandes_cv_synthese.csv',row.names=FALSE)
grid<-expand.grid(x=gx,y=gy);grid$difference<-as.vector(b$field);grid$n_eff<-rep(b$neff,length(gy));grid$difference[grid$n_eff<20]<-NA_real_;write.csv(grid,'application/champ_cv.csv',row.names=FALSE)
