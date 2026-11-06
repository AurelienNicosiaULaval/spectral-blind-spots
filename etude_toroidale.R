library(stats)
library(utils)
source('diagnostics_core.R')
rad <- function(x) sqrt(4+sin(x)^2)
# I0(r) exp(cos(x)) is the unnormalised X marginal.
draw_x_standard <- function(n) {
 result <- numeric()
 while(length(result)<n) {
  u <- runif(n,0,2*pi)
  # I0(r) <= exp(r), and cos(x) <= 1.
  prob <- exp(cos(u)-1+log(besselI(rad(u),0,expon.scaled=TRUE))+rad(u)-sqrt(5))
  result <- c(result,u[runif(n)<prob])
 }
 result[seq_len(n)]
}
draw_y <- function(x,lambda) {
 n<-length(x); out<-numeric(n); rest<-seq_len(n)
 while(length(rest)) {
  y<-runif(length(rest),0,2*pi);r<-rad(x[rest])
  keep<-runif(length(rest))<exp(2*cos(y)+lambda*sin(x[rest])*sin(y)-r)
  out[rest[keep]]<-y[keep];rest<-rest[!keep]
 }
 out
}
field_one <- function(x,y,z) {
 gx<-2*pi*(0:15)/16; gt<-2*pi*(0:31)/32
 K<-wrapped_kernel(gx,x,.7);w<-K/rowSums(K); active<-which(1/rowSums(w^2)>=20)
 A<-wrapped_kernel(gt,z,.3);B<-wrapped_kernel(gt,y,.3)
 A<-sweep(A,2,colSums(A)*2*pi/length(gt),'/');B<-sweep(B,2,colSums(B)*2*pi/length(gt),'/')
 if(!length(active)) return(matrix(0,1,length(x)))
 do.call(rbind,lapply(active,function(i) sweep(A-B,2,w[i,],'*')))
}
rows <- lapply(seq_len(100),function(i) {
 seed<-2026120000L+i;set.seed(seed);x<-draw_x_standard(240);y<-draw_y(x,1);z<-draw_y(x,-1)
 signs<-matrix(sample(c(-1,1),240*199,TRUE),240,199)
 angle<-as.matrix(expand.grid(a=c(0,pi/2),b=c(0,pi/2)))
 p<-c(canonical=sign_test(projection_features(x,y,z,angle),signs)$p,
 fourier8=sign_test(fourier_features(x,y,z,8L),signs)$p,
 field=sign_test(field_one(x,y,z),signs)$p,mmd=mmd_test(x,y,z,signs)$p)
 data.frame(rep=i,seed=seed,method=names(p),p=unname(p),reject=unname(p)<=.05)
})
out<-do.call(rbind,rows);write.csv(out,'complement/toroidal_standard.csv',row.names=FALSE)
for_grid <- function(n) {
 a<-2*pi*(0:(n-1))/n;X<-matrix(rep(a,n),n,n);Y<-matrix(rep(a,each=n),n,n)
 norm<-mean(exp(cos(a))*besselI(rad(a),0));fx<-exp(cos(X))*besselI(rad(X),0)/(2*pi*norm)
 q1<-exp(2*cos(Y)+sin(X)*sin(Y))/(2*pi*besselI(rad(X),0))
 q2<-exp(2*cos(Y)-sin(X)*sin(Y))/(2*pi*besselI(rad(X),0))
 g<-fx*(q1-q2);den<-mean(g^2)
 coeff<-vapply(seq_len(64),function(j)mean(g*sin(2*j*(X-Y))),0.0)
 c(integral_p=mean(fx*q1)*(2*pi)^2,ratio=2*sum(coeff^2)/den,energy=den)
}
a<-for_grid(256);b<-for_grid(512);stopifnot(abs(a['integral_p']-1)<1e-12,abs(b['integral_p']-1)<1e-12)
write.csv(data.frame(grid=c(256,512),integral=c(a[1],b[1]),ratio=c(a[2],b[2]),energy=c(a[3],b[3])),'complement/toroidal_energie.csv',row.names=FALSE)
