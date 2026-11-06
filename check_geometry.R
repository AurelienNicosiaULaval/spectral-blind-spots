# Contrôles indépendants des symétrisations, avec une vue visible comme contrôle positif.
library(stats)
library(utils)

script <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1])
root <- dirname(normalizePath(script))
output <- file.path(root, "resultats_geometry")
dir.create(output, showWarnings = FALSE)
grid <- 2 * pi * (0:255) / 256
reflections <- list(cos = function(u) -u, sin = function(u) pi - u)
symmetrize <- function(f, a, b, first, second) {
  u <- grid; v <- grid
  r_u <- reflections[[first]](u); r_v <- reflections[[second]](v)
  (f(u-a,v-b) + f(r_u-a,v-b) + f(u-a,r_v-b) + f(r_u-a,r_v-b)) / 4
}
set.seed(2026100911)
origins <- c(0, pi/8, pi/3, pi/4, pi/sqrt(2), runif(12,0,2*pi))
checks <- list()
for (r in 1:3) {
  k <- 2*r
  perturbation <- function(u,v) sin(k*outer(u,v,"-"))
  for (a in origins) for (first in names(reflections)) for (second in names(reflections)) {
    residual <- max(abs(symmetrize(perturbation,a,a,first,second)))
    checks[[length(checks)+1]] <- data.frame(k,origin=a,first,second,residual)
  }
  positive <- symmetrize(perturbation,pi/(2*k),0,"cos","cos")
  expected <- -outer(cos(k*grid),cos(k*grid))
  stopifnot(max(abs(positive-expected))<1e-11, max(abs(positive))>0.999)
}
checks <- do.call(rbind,checks)
stopifnot(max(checks$residual)<1e-11)
write.csv(checks,file.path(output,"common_rotations.csv"),row.names=FALSE)

angles <- c(canonical=pi/2, quarter_turn=pi/4, band6=pi/7, irrational=pi/sqrt(2))
conditioning <- do.call(rbind,lapply(names(angles),function(name) {
  do.call(rbind,lapply(1:32,function(k) {
    a <- angles[[name]]
    B <- rbind(c(1,0),c(cos(k*a),sin(k*a)))
    singular <- svd(B,nu=0,nv=0)$d
    analytic <- sqrt(max(0,1-abs(cos(k*a))))
    stopifnot(abs(min(singular)-analytic)<2e-8)
    data.frame(angle=name,frequency=k,angle_radians=a,
               minimum_singular_value=min(singular),
               analytic_minimum=analytic,determinant=det(B))
  }))
}))
write.csv(conditioning,file.path(output,"conditioning.csv"),row.names=FALSE)

# Quadrature indépendante de la variation totale et de la perte d'arc.
eta <- 0.8
conditional_tv <- integrate(function(t) abs(eta*sin(2*(0.23-t)))/(4*pi),
                            0,2*pi,subdivisions=1000L,abs.tol=1e-10,rel.tol=1e-10)$value
arc <- function(a,b) abs(atan2(sin(a-b),cos(a-b)))
risks <- vapply(c(0,0.4,pi,5),function(z) {
  integrate(function(t) arc(t,z)*(1+eta*sin(2*(0.23-t)))/(2*pi),
            0,2*pi,subdivisions=1000L,abs.tol=1e-10,rel.tol=1e-10)$value
},numeric(1))
stopifnot(abs(conditional_tv-eta/pi)<1e-8,max(abs(risks-pi/2))<1e-8)
writeLines(c(paste("Maximum common-origin residual:",max(checks$residual)),
             "Three independent-origin controls have nonzero amplitude and match exact formula.",
             paste("Conditional TV quadrature error:",abs(conditional_tv-eta/pi)),
             paste("Arc-risk quadrature error:",max(abs(risks-pi/2))),
             "Singular-value formula checked for 32 orders and four angle choices.",
             "Numerical checks corroborate stated proofs; they do not certify originality."),
           file.path(output,"checks.txt"))
cat(readLines(file.path(output,"checks.txt")),sep="\n")
