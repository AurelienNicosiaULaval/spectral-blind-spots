library(stats)
library(utils)
library(ggplot2)
source('diagnostics_core.R')
data <- read.table(gzfile('donnees/wind.txt.gz'),header=TRUE)
stopifnot(identical(names(data),c('date','hour','wind.dir')),all(is.finite(data$wind.dir)))
data <- data[order(data$date,data$hour),]
days <- sort(unique(data$date)); stopifnot(length(days)==73,all(table(data$date)==24))
rows <- split(data,data$date)
pairs <- do.call(rbind,lapply(rows,function(day) {
 stopifnot(identical(as.character(day$hour),sprintf('%02d:00',0:23)))
 data.frame(date=day$date[1],hour=seq(0,22,2),x=day$wind.dir[seq(1,23,2)],y=day$wind.dir[seq(2,24,2)])
}))
train <- pairs[pairs$date %in% days[1:50],]; test <- pairs[pairs$date %in% days[51:73],]
stopifnot(nrow(train)==600,nrow(test)==276,!any(train$date %in% test$date))
nll <- function(par) {
 k <- exp(par[1]); l <- par[2]; m <- par[3]
 r <- sqrt(k^2+(l*sin(train$x))^2)
 -sum(k*cos(train$y-m)+l*sin(train$x)*sin(train$y-m)-log(2*pi)-log(besselI(r,0,expon.scaled=TRUE))-r)
}
fits <- lapply(c(-2,0,2),function(m) optim(c(0,0,m),nll,method='L-BFGS-B',lower=c(-5,-8,-pi),upper=c(log(8),8,pi),control=list(maxit=3000)))
fits <- fits[vapply(fits,function(f) f$convergence==0 && is.finite(f$value),logical(1))]
if(!length(fits)) stop('No optimization converged')
fit <- fits[[which.min(vapply(fits,`[[`,0.0,'value'))]]
k <- exp(fit$par[1]); l <- fit$par[2]; mu <- fit$par[3]
draw <- function(x) {
 out <- numeric(length(x)); pending <- seq_along(x)
 while(length(pending)) {
  u <- runif(length(pending),0,2*pi); r <- sqrt(k^2+(l*sin(x[pending]))^2)
  accept <- runif(length(pending))<exp(k*cos(u-mu)+l*sin(x[pending])*sin(u-mu)-r)
  out[pending[accept]] <- u[accept]; pending <- pending[!accept]
 }
 out
}
set.seed(20261030); z <- replicate(50,draw(test$x)); gx <- 2*pi*(0:31)/32; gy <- 2*pi*(0:63)/64
K <- wrapped_kernel(gx,test$x,.70); w <- K/rowSums(K); n_eff <- 1/rowSums(w^2)
Ky <- wrapped_kernel(gy,test$y,.30); Ky <- sweep(Ky,2,colSums(Ky)*2*pi/length(gy),'/')
calc <- function(z) { Kz <- wrapped_kernel(gy,z,.30); Kz <- sweep(Kz,2,colSums(Kz)*2*pi/length(gy),'/'); w %*% t(Kz-Ky) }
fields <- lapply(seq_len(ncol(z)),function(j) calc(z[,j])); avg <- Reduce('+',fields)/length(fields)
grid <- expand.grid(x=gx,y=gy); grid$difference <- as.vector(avg); grid$n_eff <- rep(n_eff,length(gy)); grid$difference[grid$n_eff<20] <- NA_real_
grid$first <- as.vector(fields[[1]]); grid$last <- as.vector(fields[[50]])
dir.create('application',showWarnings=FALSE)
write.csv(pairs,'application/paires.csv',row.names=FALSE); write.csv(grid,'application/champ.csv',row.names=FALSE)
write.csv(as.data.frame(z),'application/tirages.csv',row.names=FALSE)
write.csv(data.frame(k=k,lambda=l,mu=mu,n_train=nrow(train),n_test=nrow(test),min_neff=min(n_eff),max_neff=max(n_eff),objective=fit$value,convergence=fit$convergence),'application/ajustement.csv',row.names=FALSE)
lim <- max(abs(grid$difference),na.rm=TRUE)
# Duplicate the zero boundary at 2*pi; clip at display time, not through
# scale limits, so raster cells at the circular seam are not discarded.
plot_grid <- rbind(grid,transform(grid[grid$x==0,],x=2*pi))
plot_grid <- rbind(plot_grid,transform(plot_grid[plot_grid$y==0,],y=2*pi))
extrema <- rbind(grid[which.max(grid$difference),c('x','y')],
                 grid[which.min(grid$difference),c('x','y')])
extrema$sign <- c('+','-')
p <- ggplot(plot_grid,aes(x,y,fill=difference))+geom_raster()+
 geom_contour(data=plot_grid,aes(x=x,y=y,z=difference),inherit.aes=FALSE,
              breaks=0,color='#333333',linewidth=0.3,na.rm=TRUE)+
 geom_label(data=extrema,aes(x=x,y=y,label=sign),inherit.aes=FALSE,
            fill='white',color='black',size=5,linewidth=0,
            label.padding=unit(0.04,'lines'))+
 scale_fill_gradient2(low='#2166ac',mid='white',high='#b2182b',limits=c(-lim,lim),na.value='grey80',name='Difference\nper radian')+coord_equal(xlim=c(0,2*pi),ylim=c(0,2*pi),expand=FALSE)+scale_x_continuous(breaks=c(0,pi,2*pi),labels=c('0','pi','2pi'))+scale_y_continuous(breaks=c(0,pi,2*pi),labels=c('0','pi','2pi'))+labs(x='Predictor direction (radians)',y='Response direction (radians)')+theme_minimal(base_size=11)
ggsave('application/champ_vent.pdf',p,width=7,height=5)
writeLines(capture.output(sessionInfo()),'application/sessionInfo.txt')
