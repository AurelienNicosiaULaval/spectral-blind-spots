library(stats)
library(utils)
p <- read.csv('application/paires.csv'); first <- sort(unique(p$date))[1:50]; d <- p[p$date %in% first,]
f <- function(v) {
 k<-exp(v[1]);lambda<-v[2];mu<-v[3]; rad<-sqrt(k*k+(lambda*sin(d$x))^2)
 -sum(k*cos(d$y-mu)+lambda*sin(d$x)*sin(d$y-mu)-log(2*pi)-log(besselI(rad,0,expon.scaled=TRUE))-rad)
}
starts <- expand.grid(logk=c(-1,0,1),lambda=c(-2,0,2),mu=c(-2,0,2))
fits <- lapply(seq_len(nrow(starts)),function(i) optim(as.numeric(starts[i,]),f,method='L-BFGS-B',lower=c(-5,-8,-pi),upper=c(log(8),8,pi),control=list(maxit=3000)))
r <- data.frame(start=seq_along(fits),value=vapply(fits,`[[`,0.0,'value'),convergence=vapply(fits,`[[`,0L,'convergence'))
write.csv(r,'application/verification_optimisation.csv',row.names=FALSE)
expected<-read.csv('application/ajustement.csv')$objective
if(abs(min(r$value[r$convergence==0])-expected)>1e-5) stop('A better fit was found; inspect before relying on example')
cat('Best over',nrow(starts),'starts matches objective',expected,'\n')
