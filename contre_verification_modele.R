library(stats)
library(utils)
# Independent grid integration of conditional moments against rejection draws.
grid <- 2*pi*((0:8191)+.5)/8192
sample_y <- function(n,x,k,lambda,mu) {
  out<-numeric(); radius<-sqrt(k^2+(lambda*sin(x))^2)
  while(length(out)<n) {
    z<-runif(n,0,2*pi);u<-runif(n)
    out<-c(out,z[u<exp(k*cos(z-mu)+lambda*sin(x)*sin(z-mu)-radius)])
  }
  out[seq_len(n)]
}
set.seed(202610301L)
rows <- lapply(c(0,pi/2,pi,3*pi/2),function(x) {
  weights<-exp(2*cos(grid)+sin(x)*sin(grid));weights<-weights/sum(weights)
  z<-sample_y(20000,x,2,1,0)
  data.frame(x=x,cos_ref=sum(weights*cos(grid)),cos_mc=mean(cos(z)),cos_se=sd(cos(z))/sqrt(length(z)),
    sin_ref=sum(weights*sin(grid)),sin_mc=mean(sin(z)),sin_se=sd(sin(z))/sqrt(length(z)))
})
out<-do.call(rbind,rows)
stopifnot(all(abs(out$cos_ref-out$cos_mc)<5*out$cos_se),all(abs(out$sin_ref-out$sin_mc)<5*out$sin_se))
write.csv(out,'complement/validation_moments.csv',row.names=FALSE)
cat('Grid-integrated moments agree with rejection draws to within five Monte Carlo standard errors at four predictors.\n')
