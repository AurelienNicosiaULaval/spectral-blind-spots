library(stats)
source('diagnostics_core.R')
draw_vm <- function(x) {
 out <- numeric(length(x)); todo <- seq_along(x)
 while(length(todo)) {
  u<-runif(length(todo),0,2*pi); r<-sqrt(4+sin(x[todo])^2)
  keep<-runif(length(todo))<exp(2*cos(u)+sin(x[todo])*sin(u)-r)
  out[todo[keep]]<-u[keep]; todo<-todo[!keep]
 }
 out
}
results <- lapply(seq_len(500),function(i) {
 seed<-2026111000L+i;set.seed(seed); x<-runif(240,0,2*pi);y<-draw_vm(x);z<-draw_vm(x)
 signs<-matrix(sample(c(-1,1),240*199,replace=TRUE),240,199)
 p<-sign_test(fourier_features(x,y,z,12L),signs)$p
 data.frame(rep=i,seed=seed,p=p,reject=p<=.05)
})
out<-do.call(rbind,results);write.csv(out,'complement/calibration_supplementaire.csv',row.names=FALSE)
