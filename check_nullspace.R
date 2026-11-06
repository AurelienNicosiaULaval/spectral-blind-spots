# Independent finite-block checks for the complete common-origin null-space proof.
library(stats)
library(utils)
script <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1])
root <- dirname(normalizePath(script))
origins <- 2*pi*(0:100)/101
rows <- list()
for (k in 1:12) for (l in 1:12) {
  # The number of common origins resolves all sums/differences through order 24.
  map <- do.call(rbind,lapply(origins,function(a) {
    do.call(rbind,lapply(0:1,function(j) {
      do.call(rbind,lapply(0:1,function(m) {
        kronecker(c(cos(k*(a+j*pi/2)),sin(k*(a+j*pi/2))),
                  c(cos(l*(a+m*pi/2)),sin(l*(a+m*pi/2))))
      }))
    }))
  }))
  singular <- svd(map,nu=0,nv=4)
  observed_rank <- sum(singular$d>1e-9)
  expected_rank <- if(k==l && k%%2==0) 3L else 4L
  stopifnot(observed_rank==expected_rank)
  null_residual <- NA_real_
  if(expected_rank==3L) {
    antisymmetric <- c(0,-1,1,0)/sqrt(2)
    null_residual <- max(abs(map%*%antisymmetric))
    stopifnot(null_residual<1e-11,
              abs(sum(singular$v[,4]*antisymmetric))>1-1e-10)
  }
  rows[[length(rows)+1L]] <- data.frame(k,l,observed_rank,expected_rank,
                                       null_residual)
}
result <- do.call(rbind,rows)
write.csv(result,file.path(root,"resultats_geometry","nullspace_blocks.csv"),row.names=FALSE)
# Finite set of four common origins used by the simulation (0, pi/8, pi/3, pi/4), orders 0 through 8.
# Added in the 10 October 2026 revision: a finite origin set can have blind directions beyond Proposition 2.
finite_origins <- c(0,pi/8,pi/3,pi/4)
direction <- function(k,a) if(k==0) 1 else c(cos(k*a),sin(k*a))
finite <- list()
for (k in 0:8) for (l in 0:8) {
  if (k==0 && l==0) next
  map <- do.call(rbind,lapply(finite_origins,function(a) {
    do.call(rbind,lapply(0:1,function(j) do.call(rbind,lapply(0:1,function(m) {
      kronecker(direction(k,a+j*pi/2),direction(l,a+m*pi/2))
    }))))
  }))
  singular <- svd(map,nu=0,nv=0)$d
  nullity <- ncol(map)-sum(singular>1e-9)
  all_origin_nullity <- as.integer(k==l && k>0 && k%%2==0)
  finite[[length(finite)+1L]] <- data.frame(k,l,dimension=ncol(map),nullity,
                                           all_origin_nullity,extra=nullity-all_origin_nullity)
}
finite <- do.call(rbind,finite)
stopifnot(sum(finite$nullity)==11L,sum(finite$extra)==7L,
          all(finite$extra[finite$extra>0] > 0 & (finite$k[finite$extra>0]==8 | finite$l[finite$extra>0]==8)))
write.csv(finite,file.path(root,"resultats_geometry","finite_origins_blocks.csv"),row.names=FALSE)
extra_blocks <- finite[finite$extra>0,]
writeLines(c("144 real Fourier tensor blocks checked, coordinate orders 1 through 12.",
             "Exactly the six equal even-order blocks have rank three; all others have rank four.",
             paste("Maximum residual for the predicted antisymmetric null vector:",
                   max(result$null_residual,na.rm=TRUE)),
             paste("Four simulation origins 0, pi/8, pi/3, pi/4, orders 0-8: blind dimension",sum(finite$nullity),
                   "= 4 from Proposition 2 +",sum(finite$extra),"additional directions in blocks",
                   paste(sprintf("(%d,%d)",extra_blocks$k,extra_blocks$l),collapse=" ")),
             "The additional directions arise because 8a is a multiple of pi for a = pi/8 and pi/4.",
             "This finite calculation corroborates the analytical L2 proof; it is not that proof or an originality certificate."),
           file.path(root,"resultats_geometry","nullspace_checks.txt"))
cat(readLines(file.path(root,"resultats_geometry","nullspace_checks.txt")),sep="\n")
