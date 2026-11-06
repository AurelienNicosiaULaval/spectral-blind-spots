library(stats)
library(utils)

wrap <- function(x) x %% (2*pi)
arc <- function(a,b) abs(atan2(sin(a-b),cos(a-b)))
bump <- function(x,center) pmax(1-(arc(x,center)/0.9)^2,0)^2
draw_x <- function(n,rho) {
  result <- numeric(0)
  while(length(result)<n) {
    proposal <- runif(n,0,2*pi)
    result <- c(result,proposal[runif(n)<(1+rho*cos(proposal))/(1+abs(rho))])
  }
  result[seq_len(n)]
}
perturbation <- function(x,theta,kind,order,amplitude) {
  if(kind=="null_uniform") return(rep(0,length(theta)))
  if(kind %in% c("phase","null_conditional")) return(amplitude*sin(order*(x-theta)))
  if(kind=="product") return(amplitude*sin(2*x)*sin(2*theta))
  if(kind=="local_rare") return(amplitude*bump(x,pi)*sin(2*theta))
  if(kind=="local_common") return(amplitude*bump(x,0)*sin(2*theta))
  stop("Unknown data-generating law")
}
draw_response <- function(x,definition) {
  result <- numeric(length(x)); pending <- seq_along(x)
  while(length(pending)) {
    proposal <- runif(length(pending),0,2*pi)
    p <- perturbation(x[pending],proposal,definition$kind,
                      definition$order,definition$amplitude)
    keep <- runif(length(pending))<(1+p)/(1+definition$amplitude)
    result[pending[keep]] <- proposal[keep]
    pending <- pending[!keep]
  }
  result
}
sign_test <- function(contributions,signs) {
  denominator <- sqrt(rowSums(contributions^2))
  active <- denominator>max(denominator)*1e-12
  if(!any(active)) return(list(p=1,statistic=0,components=0L))
  standardized <- contributions[active,,drop=FALSE]/denominator[active]
  observed <- max(abs(rowSums(standardized)))
  randomized <- apply(abs(standardized %*% signs),2,max)
  p <- (1+sum(randomized>=observed-1e-12*max(1,observed)))/(ncol(signs)+1)
  list(p=p,statistic=observed,components=sum(active))
}
projection_features <- function(x,y,z,angle_pairs,degree=8L) {
  do.call(rbind,lapply(seq_len(nrow(angle_pairs)),function(index) {
    a <- angle_pairs[index,1]; b <- angle_pairs[index,2]
    X <- cos(outer(0:degree,x-a))
    difference <- cos(outer(1:degree,z-b))-cos(outer(1:degree,y-b))
    do.call(rbind,lapply(seq_len(nrow(X)),function(k) {
      sweep(difference,2,X[k,],"*")
    }))
  }))
}
fourier_features <- function(x,y,z,degree=8L) {
  frequencies <- expand.grid(k=-degree:degree,l=1:degree)
  values <- vapply(seq_len(nrow(frequencies)),function(index) {
    exp(1i*(frequencies$k[index]*x+frequencies$l[index]*z))-
      exp(1i*(frequencies$k[index]*x+frequencies$l[index]*y))
  },complex(length(x)))
  rbind(t(Re(values)),t(Im(values)))
}
wrapped_kernel <- function(a,b,h) {
  d <- outer(wrap(a),wrap(b),"-")
  Reduce(`+`,lapply(-2:2,function(index) dnorm(d+2*pi*index,sd=h)))
}
field_features <- function(x,y,z) {
  gx <- 2*pi*(0:15)/16; gt <- 2*pi*(0:31)/32
  blocks <- list(); min_support <- Inf; max_support <- 0L; identity_error <- 0
  for(hx in c(0.35,0.70)) {
    K <- wrapped_kernel(gx,x,hx); weights <- K/rowSums(K)
    effective <- 1/rowSums(weights^2); selected <- which(effective>=20)
    min_support <- min(min_support,length(selected)); max_support <- max(max_support,length(selected))
    for(ht in c(0.15,0.30)) {
      Ky <- wrapped_kernel(gt,y,ht); Kz <- wrapped_kernel(gt,z,ht)
      Ky <- sweep(Ky,2,colSums(Ky)*2*pi/length(gt),"/")
      Kz <- sweep(Kz,2,colSums(Kz)*2*pi/length(gt),"/")
      difference <- Kz-Ky
      conditional <- weights %*% t(difference)
      witness <- (K/length(x)) %*% t(difference)
      identity_error <- max(identity_error,max(abs(conditional-witness/rowMeans(K))))
      if(length(selected)) blocks[[length(blocks)+1]] <- do.call(rbind,lapply(selected,function(i) {
        sweep(difference,2,weights[i,],"*")
      }))
    }
  }
  if(!length(blocks)) return(list(C=matrix(0,1,length(x)),min_support=0L,
                                max_support=0L,identity_error=identity_error))
  list(C=do.call(rbind,blocks),min_support=min_support,max_support=max_support,
       identity_error=identity_error)
}
mmd_test <- function(x,y,z,signs) {
  kernel <- function(a,b) exp(-(1-cos(outer(a,b,"-")))/0.5^2)
  Ky <- kernel(y,y); Kz <- kernel(z,z); Kyz <- kernel(y,z)
  H <- kernel(x,x)*(Ky+Kz-Kyz-t(Kyz))
  observed <- sum(H)
  randomized <- colSums(signs*(H %*% signs))
  p <- (1+sum(randomized>=observed-1e-11*max(1,abs(observed))))/(ncol(signs)+1)
  list(p=p,statistic=observed/length(x)^2,components=length(x)^2)
}
candidate_pit <- function(theta,x,definition) {
  if(!definition$is_null || definition$kind=="null_uniform") return(theta/(2*pi))
  # Intégrale de sin(x-t) entre 0 et theta, pour la nullité conditionnelle d'ordre 1.
  (theta+definition$amplitude*(cos(x-theta)-cos(x)))/(2*pi)
}
coverage_test <- function(x,y,definition,B=199L) {
  u <- candidate_pit(y,x,definition)
  probabilities <- seq(0.05,0.95,by=0.05)
  gx <- 2*pi*(0:15)/16
  null_u <- matrix(runif(length(x)*B),length(x),B)
  observed <- numeric(0); randomized <- matrix(numeric(0),0,B)
  for(hx in c(0.35,0.70)) {
    K <- wrapped_kernel(gx,x,hx); w <- K/rowSums(K)
    included <- which(1/rowSums(w^2)>=20)
    if(!length(included)) next
    w <- w[included,,drop=FALSE]
    for(alpha in probabilities) {
      sd <- sqrt(alpha*(1-alpha)*rowSums(w^2))
      observed <- c(observed,as.vector(w %*% (u<=alpha)-alpha)/sd)
      randomized <- rbind(randomized,(w %*% (null_u<=alpha)-alpha)/sd)
    }
  }
  if(!length(observed)) return(list(p=1,statistic=0,components=0L))
  statistic <- max(abs(observed)); reference <- apply(abs(randomized),2,max)
  list(p=(1+sum(reference>=statistic-1e-12*max(1,statistic)))/(B+1),
       statistic=statistic,components=length(observed))
}
