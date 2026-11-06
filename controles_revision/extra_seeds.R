# Additional datasets with fresh seeds (outside the study's seed range), using the authors' core functions,
# to check the population-blind cells of the canonical/common/crossed_pi4 diagnostics at n=240 (phase 4).
args <- commandArgs(TRUE); core <- args[1]; out <- args[2]
source(core)
angles <- list(canonical=as.matrix(expand.grid(a=c(0,pi/2),b=c(0,pi/2))),
  crossed_pi4=as.matrix(expand.grid(a=c(0,pi/4),b=c(0,pi/4))))
angles$common_rotations <- do.call(rbind,lapply(c(0,pi/8,pi/3,pi/4),function(a) angles$canonical+a))
res <- list()
for (cfg in list(list(order=4,n=240L),list(order=6,n=240L))) {
  def <- list(kind="phase",order=cfg$order,amplitude=0.8,is_null=FALSE)
  for (r in 1:1000) {
    set.seed(1500000000L + 100000L*as.integer(cfg$order) + r)
    x <- draw_x(cfg$n,0); y <- draw_response(x,def); z <- runif(cfg$n,0,2*pi)
    signs <- matrix(sample(c(-1,1),cfg$n*199,TRUE),cfg$n,199)
    for (nm in names(angles)) res[[length(res)+1]] <- data.frame(order=cfg$order,n=cfg$n,rep=r,diagnostic=nm,
       p=sign_test(projection_features(x,y,z,angles[[nm]]),signs)$p)
  }
}
res <- do.call(rbind,res); write.csv(res,out,row.names=FALSE)
print(aggregate(p~order+diagnostic,res,function(p) c(rate=mean(p<=0.05),k=sum(p<=0.05))))
