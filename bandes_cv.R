library(stats)
library(utils)
source('diagnostics_core.R')
p<-read.csv('application/paires.csv');days<-sort(unique(p$date));train<-p[p$date %in% days[1:50],]
z<-as.matrix(read.csv('application/tirages.csv'));test<-p[p$date %in% days[51:73],]
# Stable log-sums retain finite scores even on the smaller training-only grid.
logadd<-function(a,b){m<-pmax(a,b);m+log(exp(a-m)+exp(b-m))}
logkernel<-function(a,b,h) {
 d<-outer(wrap(a),wrap(b),'-')
 Reduce(logadd,lapply(-2:2,function(j)dnorm(d+2*pi*j,sd=h,log=TRUE)))
}
rowlse<-function(M){m<-apply(M,1,max);m+log(rowSums(exp(M-m)))}
hx<-c(.0875,.175,.35,.7,1.2);ht<-c(.0375,.075,.15,.3,.6)
LX<-lapply(hx,function(h)logkernel(train$x,train$x,h));LY<-lapply(ht,function(h)logkernel(train$y,train$y,h))
excluded<-outer(train$date,train$date,'==');LX<-lapply(LX,function(m){m[excluded]<--Inf;m})
cv<-expand.grid(hx=hx,ht=ht)
cv$loglik<-vapply(1:nrow(cv),function(i){
 a<-LX[[match(cv$hx[i],hx)]];b<-LY[[match(cv$ht[i],ht)]]
 score<-rowlse(a+b)-rowlse(a);stopifnot(all(is.finite(score)));sum(score)
},numeric(1))
cv$selected<-seq_len(nrow(cv))==which.max(cv$loglik)
cv$boundary<-cv$hx%in%range(hx)|cv$ht%in%range(ht)
write.csv(cv,'application/bandes_cv.csv',row.names=FALSE)
best<-cv[cv$selected,]
gx<-2*pi*(0:31)/32;gy<-2*pi*(0:63)/64
field<-function(hx,ht){
 K<-wrapped_kernel(gx,test$x,hx);w<-K/rowSums(K);neff<-1/rowSums(w^2)
 ky<-wrapped_kernel(gy,test$y,ht);ky<-sweep(ky,2,colSums(ky)*2*pi/length(gy),'/')
 f<-matrix(0,length(gx),length(gy))
 for(j in seq_len(ncol(z))){kz<-wrapped_kernel(gy,z[,j],ht);kz<-sweep(kz,2,colSums(kz)*2*pi/length(gy),'/');f<-f+w%*%t(kz-ky)/ncol(z)}
 list(field=f,neff=neff)
}
a<-field(.7,.3);b<-field(best$hx,best$ht);mask<-a$neff>=20&b$neff>=20
stopifnot(any(mask))
summary<-data.frame(hx=best$hx,ht=best$ht,loglik=best$loglik,boundary=best$boundary,active=sum(mask),minimum_support=min(b$neff),rms_diff=sqrt(mean((a$field[mask,]-b$field[mask,])^2)),max_abs_diff=max(abs(a$field[mask,]-b$field[mask,])))
write.csv(summary,'application/bandes_cv_synthese.csv',row.names=FALSE)
grid<-expand.grid(x=gx,y=gy);grid$difference<-as.vector(b$field);grid$n_eff<-rep(b$neff,length(gy));grid$difference[grid$n_eff<20]<-NA_real_
write.csv(grid,'application/champ_cv.csv',row.names=FALSE)
print(best);print(summary)
