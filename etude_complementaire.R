library(stats)
library(utils)
source('diagnostics_core.R')

# New, post-review study. The reference and candidate share a uniform X.
cases <- data.frame(case=c('fidelity','direction','concentration','association','harmonic9'),
  kr=c(2,2,2,2,0), lr=c(1,1,1,1,0), mr=c(0,0,0,0,0),
  kq=c(2,2,1,2,0), lq=c(1,1,1,-1,0), mq=c(0,.5,0,0,0))
draw_conditional <- function(x,k,l,m) {
  n <- length(x); y <- numeric(n); todo <- seq_len(n)
  while(length(todo)) {
    u <- runif(length(todo),0,2*pi)
    radius <- sqrt(k^2+(l*sin(x[todo]))^2)
    accept <- runif(length(todo)) < exp(k*cos(u-m)+l*sin(x[todo])*sin(u-m)-radius)
    y[todo[accept]] <- u[accept]; todo <- todo[!accept]
  }
  y
}
density_joint <- function(x,y,k,l,m) {
  radius <- sqrt(k^2+(l*sin(x))^2)
  exp(k*cos(y-m)+l*sin(x)*sin(y-m))/(4*pi^2*besselI(radius,0))
}
blind_energy <- function(case, grid=256L, maxorder=64L) {
  x <- (0:(grid-1))*2*pi/grid; theta <- x
  X <- matrix(rep(x,grid),grid,grid); Y <- matrix(rep(theta,each=grid),grid,grid)
  if(case$case=='harmonic9') {
    g <- .8*sin(9*(X-Y))/(4*pi^2)
  } else {
    g <- density_joint(X,Y,case$kr,case$lr,case$mr)-density_joint(X,Y,case$kq,case$lq,case$mq)
  }
  denom <- mean(g^2)
  if(denom < 1e-30) return(c(energy=denom,blind=NA_real_))
  coeff <- vapply(seq_len(maxorder),function(j) mean(g*sin(2*j*(X-Y))),0.0)
  c(energy=denom,blind=2*sum(coeff^2)/denom)
}
field_single <- function(x,y,z,hx,ht) {
  gx <- 2*pi*(0:15)/16; gt <- 2*pi*(0:31)/32
  K <- wrapped_kernel(gx,x,hx); w <- K/rowSums(K)
  active <- which(1/rowSums(w^2)>=20)
  if(!length(active)) return(matrix(0,1,length(x)))
  A <- wrapped_kernel(gt,z,ht); B <- wrapped_kernel(gt,y,ht)
  A <- sweep(A,2,colSums(A)*2*pi/length(gt),'/')
  B <- sweep(B,2,colSums(B)*2*pi/length(gt),'/')
  do.call(rbind,lapply(active,function(i) sweep(A-B,2,w[i,],'*')))
}
mmd_scale <- function(x,y,z,signs,scale) {
  kern <- function(a,b) exp(-(1-cos(outer(a,b,'-')))/scale^2)
  H <- kern(x,x)*(kern(y,y)+kern(z,z)-kern(y,z)-t(kern(y,z)))
  obs <- sum(H); sim <- colSums(signs*(H %*% signs))
  (1+sum(sim>=obs-1e-11*max(1,abs(obs))))/(ncol(signs)+1)
}
rows <- list(); failures <- list(); ix <- 0L; set.seed(20261030L)
for(i in seq_len(nrow(cases))) for(rep in seq_len(100L)) {
  scenario <- cases[i,]; seed <- 2026103001L+1000L*i+rep
  result <- tryCatch({
    set.seed(seed); n <- 240L; x <- runif(n,0,2*pi)
    if(scenario$case=='harmonic9') {
      def <- list(kind='phase',order=9,amplitude=.8)
      y <- draw_response(x,def); z <- runif(n,0,2*pi)
    } else {
      y <- draw_conditional(x,scenario$kr,scenario$lr,scenario$mr)
      z <- draw_conditional(x,scenario$kq,scenario$lq,scenario$mq)
    }
    signs <- matrix(sample(c(-1,1),n*199,replace=TRUE),n,199)
    angles <- as.matrix(expand.grid(a=c(0,pi/2),b=c(0,pi/2)))
    p <- c(canonical=sign_test(projection_features(x,y,z,angles),signs)$p,
      fourier4=sign_test(fourier_features(x,y,z,4L),signs)$p,
      fourier8=sign_test(fourier_features(x,y,z,8L),signs)$p,
      fourier12=sign_test(fourier_features(x,y,z,12L),signs)$p,
      field_narrow=sign_test(field_single(x,y,z,.35,.15),signs)$p,
      field_wide=sign_test(field_single(x,y,z,.70,.30),signs)$p,
      mmd03=mmd_scale(x,y,z,signs,.3),mmd05=mmd_scale(x,y,z,signs,.5),
      mmd10=mmd_scale(x,y,z,signs,1))
    data.frame(case=scenario$case,rep=rep,seed=seed,method=names(p),p=unname(p),reject=unname(p)<=.05)
  },error=function(e) {failures[[length(failures)+1L]] <<- data.frame(case=scenario$case,rep=rep,seed=seed,error=conditionMessage(e)); NULL})
  if(!is.null(result)) {ix <- ix+1L; rows[[ix]] <- result}
}
dir.create('complement',showWarnings=FALSE)
write.csv(if(length(rows)) do.call(rbind,rows) else data.frame(), 'complement/resultats.csv',row.names=FALSE)
write.csv(if(length(failures)) do.call(rbind,failures) else data.frame(case=character(),rep=integer(),seed=integer(),error=character()),'complement/echecs.csv',row.names=FALSE)
if(length(failures)) stop('Failures recorded; inspect before reporting.')
raw <- do.call(rbind,rows)
groups <- split(raw,interaction(raw$case,raw$method,drop=TRUE))
summary <- do.call(rbind,lapply(groups,function(d) {
  a <- sum(d$reject); n <- nrow(d); z <- qnorm(.975); den <- 1+z^2/n
  mid <- (a/n+z^2/(2*n))/den; hw <- z*sqrt((a/n)*(1-a/n)/n+z^2/(4*n^2))/den
  data.frame(case=d$case[1],method=d$method[1],n=n,rejections=a,rate=a/n,lower=mid-hw,upper=mid+hw)
}))
write.csv(summary,'complement/synthese.csv',row.names=FALSE)
energy <- do.call(rbind,lapply(seq_len(nrow(cases)),function(i) {
  a <- blind_energy(cases[i,],256L); b <- blind_energy(cases[i,],512L)
  data.frame(case=cases$case[i],energy=a['energy'],ratio256=a['blind'],ratio512=b['blind'])
}))
write.csv(energy,'complement/energie.csv',row.names=FALSE)
writeLines(capture.output(sessionInfo()),'complement/sessionInfo.txt')
