library(stats)
library(utils)
source('diagnostics_core.R')
# For each row and each permutation choose which of m+1 responses is the
# reference; recompute all feature normalizers for the new labelling.
multiple_test <- function(x,y,z,permutation,B) {
 n<-length(x);m<-ncol(z);u<-cbind(y,z);deg<-4L
 freq<-expand.grid(k=(-deg):deg,l=1:deg)
 phase<-outer(freq$k,x,'*')
 tensor<-array(0+0i,c(nrow(freq),n,m+1L))
 for(j in seq_len(m+1L)) tensor[,,j]<-exp(1i*(phase+outer(freq$l,u[,j],'*')))
 totals<-apply(tensor,c(1,2),sum)
 make <- function(idx) {
  chosen<-matrix(tensor[cbind(rep(seq_len(nrow(freq)),times=n),rep(seq_len(n),each=nrow(freq)),rep(idx,each=nrow(freq)))],nrow=nrow(freq))
  v<-(totals-(m+1)*chosen)/m
  rbind(Re(v),Im(v))
 }
 statistic <- function(C) {
  d<-sqrt(rowSums(C^2));active<-d>max(d)*1e-12
  if(!any(active)) return(0)
  max(abs(rowSums(C[active,,drop=FALSE])/d[active]))
 }
 obs<-statistic(make(rep(1L,n)));random<-vapply(seq_len(B),function(b)statistic(make(permutation[,b])),0.0)
 (1+sum(random>=obs-1e-12*max(1,obs)))/(B+1)
}
scenarios<-c('fidelity','phase2');rows<-list(); failures<-list();ix<-0L
for(i in seq_along(scenarios)) for(rep in seq_len(100)) {
 seed<-2026100010L+i*1000L+rep
 ans<-tryCatch({
  set.seed(seed);n<-240L;x<-runif(n,0,2*pi)
  def<-if(i==1) list(kind='null_uniform',order=0,amplitude=0,is_null=TRUE) else list(kind='phase',order=2,amplitude=.4,is_null=FALSE)
  y<-draw_response(x,def); z<-replicate(5,runif(n,0,2*pi)); B<-199L
  signs<-matrix(sample(c(-1,1),n*B,TRUE),n,B)
  one<-sign_test(fourier_features(x,y,z[,1],4L),signs)$p
  labels<-matrix(sample.int(6L,n*B,TRUE),n,B)
  five<-multiple_test(x,y,z,labels,B)
  # Each transformed angle is a CDF at the given cut for the uniform candidate.
  pit<-vapply(c(0,pi/2,pi),function(a) coverage_test(x,wrap(y-a),def,B)$p,0.0)
  p<-c(one=one,five=five,cut0=pit[1],cutpi2=pit[2],cutpi=pit[3])
  data.frame(scenario=scenarios[i],rep=rep,seed=seed,method=names(p),p=unname(p),reject=unname(p)<=.05)
 },error=function(e){ failures[[length(failures)+1]]<<-data.frame(scenario=scenarios[i],rep=rep,error=conditionMessage(e));NULL})
 if(!is.null(ans)){ix<-ix+1L;rows[[ix]]<-ans}
}
dir.create('complement',showWarnings=FALSE)
write.csv(if(length(failures)) do.call(rbind,failures) else data.frame(scenario=character(),rep=integer(),error=character()),'complement/echecs_origine_tirages.csv',row.names=FALSE)
if(length(failures)) stop('Failures recorded')
out<-do.call(rbind,rows);write.csv(out,'complement/origine_tirages.csv',row.names=FALSE)
