library(stats)
library(utils)
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")][1])
root <- dirname(normalizePath(script))
source(file.path(root,"diagnostics_core.R"))

# Contrôle de la calibration sur un groupe complet de 256 échanges, indépendant du générateur.
set.seed(481)
C <- matrix(rnorm(5*8),5,8)
signs <- t(as.matrix(expand.grid(rep(list(c(-1,1)),8))))
p <- vapply(seq_len(ncol(signs)),function(i) {
  sign_test(sweep(C,2,signs[,i],"*"),signs)$p
},numeric(1))
stopifnot(mean(p<=0.05)<=0.05,mean(p<=0.10)<=0.10)

# Le champ et la matrice quadratique se transforment correctement sous un échange.
x <- runif(32,0,2*pi); y <- runif(32,0,2*pi); z <- runif(32,0,2*pi)
epsilon <- sample(c(-1,1),32,TRUE)
y_swap <- ifelse(epsilon==1,y,z); z_swap <- ifelse(epsilon==1,z,y)
field <- field_features(x,y,z); swapped <- field_features(x,y_swap,z_swap)
stopifnot(field$identity_error<1e-12,
          max(abs(swapped$C-sweep(field$C,2,epsilon,"*")))<1e-12)
kernel <- function(a,b) exp(-(1-cos(outer(a,b,"-")))/0.5^2)
H <- kernel(x,x)*(kernel(y,y)+kernel(z,z)-kernel(y,z)-kernel(z,y))
H_swap <- kernel(x,x)*(kernel(y_swap,y_swap)+kernel(z_swap,z_swap)-
                       kernel(y_swap,z_swap)-kernel(z_swap,y_swap))
stopifnot(max(abs(H_swap-H*outer(epsilon,epsilon)))<1e-12)

# Validation de la CDF utilisée par le diagnostic supplémentaire, avec intégration indépendante.
definition <- list(is_null=TRUE,kind="null_conditional",amplitude=0.8,order=1)
errors <- vapply(c(0.2,1.4,4.7),function(x0) {
  density <- function(theta) (1+0.8*sin(x0-theta))/(2*pi)
  mean_u <- integrate(function(theta) candidate_pit(theta,x0,definition)*density(theta),
                      0,2*pi,abs.tol=1e-10,rel.tol=1e-10)$value
  second_u <- integrate(function(theta) candidate_pit(theta,x0,definition)^2*density(theta),
                        0,2*pi,abs.tol=1e-10,rel.tol=1e-10)$value
  max(abs(mean_u-0.5),abs(second_u-1/3))
},numeric(1))
stopifnot(max(errors)<1e-9)
writeLines(c("Full exchange-group calibration checked at 5% and 10%.",
             paste("Field exchange-equivariance residual:",
                   max(abs(swapped$C-sweep(field$C,2,epsilon,"*")))),
             paste("Quadratic exchange-equivariance residual:",
                   max(abs(H_swap-H*outer(epsilon,epsilon)))),
             paste("Candidate CDF moment quadrature residual:",max(errors))),
           file.path(root,"checks_diagnostics.txt"))
cat(readLines(file.path(root,"checks_diagnostics.txt")),sep="\n")
