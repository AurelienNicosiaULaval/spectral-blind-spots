library(stats)
library(utils)
dir.create('renforcement',showWarnings=FALSE)
g_wc<-function(t,rho)(1-rho^2)/(2*pi*(1+rho^2-2*rho*cos(t)))
energy_grid<-function(rho,delta,n,axial=FALSE) {
 a<-2*pi*(0:(n-1))/n;X<-matrix(rep(a,n),n,n);Y<-matrix(rep(a,each=n),n,n)
 f<-if(axial)2 else 1
 pp<-g_wc(f*(Y-X-delta),rho)/(2*pi);qq<-g_wc(f*(Y-X+delta),rho)/(2*pi)
 d<-pp-qq
 coef<-vapply(1:64,function(j)mean(d*sin(2*j*(X-Y))),numeric(1))
 c(ratio=2*sum(coef^2)/mean(d^2),integral=mean(pp)*(2*pi)^2,
   margin_error=max(abs(rowMeans(pp)*2*pi-1/(2*pi)),abs(colMeans(pp)*2*pi-1/(2*pi))))
}
cases<-expand.grid(rho=c(.5,.8,.9),delta=c(.2,.5,1));cases$axial<-FALSE
cases<-rbind(cases,data.frame(rho=.8,delta=.5,axial=TRUE))
rows<-lapply(1:nrow(cases),function(i) {
 ca<-cases[i,];k<-1:500;f<-if(ca$axial)2 else 1;w<-ca$rho^(2*k)*sin(f*k*ca$delta)^2
 exact<-if(ca$axial)1 else sum(w[k%%2==0])/sum(w)
 a<-energy_grid(ca$rho,ca$delta,256,ca$axial);b<-energy_grid(ca$rho,ca$delta,512,ca$axial)
 stopifnot(abs(a['ratio']-exact)<1e-8,abs(b['ratio']-exact)<1e-8,abs(a['integral']-1)<1e-10,abs(b['integral']-1)<1e-10,a['margin_error']<1e-10,b['margin_error']<1e-10)
 cbind(ca,data.frame(series=exact,grid256=a['ratio'],grid512=b['ratio'],integral=b['integral'],margin_error=b['margin_error']))
})
write.csv(do.call(rbind,rows),'renforcement/rotation_energie.csv',row.names=FALSE)
cat('Series, quadrature, unit integrals and uniform margins verified for all 10 cases.\n')
