library(stats)
library(utils)
library(grDevices)
library(grid)
library(ggplot2)
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")][1])
root <- dirname(normalizePath(script))
out <- file.path(root,"figures");dir.create(out,showWarnings=FALSE)
eta <- 0.8
angular <- expand.grid(x=2*pi*(0:95)/96,theta=2*pi*(0:95)/96)
angular$density <- (1+eta*sin(2*(angular$x-angular$theta)))/(2*pi)
reference <- ggplot(angular,aes(x,theta,fill=density))+
  geom_tile(width=2*pi/96,height=2*pi/96)+
  scale_fill_gradient(low="#FAFAF7",high="#27649A",name="Density (1/radian)",
    guide=guide_colourbar(title.position="top",title.hjust=0.5,barwidth=unit(2.2,"cm")))+
  scale_x_continuous(breaks=c(0,pi,2*pi),labels=c("0","pi","2pi"),expand=expansion(mult=0))+
  scale_y_continuous(breaks=c(0,pi,2*pi),labels=c("0","pi","2pi"),expand=expansion(mult=0))+
  coord_cartesian(xlim=c(0,2*pi),ylim=c(0,2*pi),expand=FALSE)+
  labs(x="Predictor (radians)",y="Response (radians)",title="Reference conditional law")+
  theme_minimal(base_size=12)+theme(panel.grid=element_blank(),legend.position="bottom",
                                  plot.title=element_text(size=12))
projected <- expand.grid(s=seq(-0.99,0.99,length.out=96),t=seq(-0.99,0.99,length.out=96))
make_projection <- function(a,b,title) {
  projected$relative_difference <- -eta*sin(2*(a-b))*(2*projected$s^2-1)*(2*projected$t^2-1)
  p <- ggplot(projected,aes(s,t,fill=relative_difference))+
    geom_raster()+scale_fill_gradient2(low="#27649A",mid="#FAFAF7",high="#AC3A32",
      midpoint=0,limits=c(-eta,eta),breaks=c(-eta,0,eta),name="Relative density difference",
      guide=guide_colourbar(title.position="top",title.hjust=0.5,barwidth=unit(2.2,"cm")))+
    coord_equal()+labs(x="Cosine of predictor",y="Cosine of response",title=title)+
    theme_minimal(base_size=12)+theme(panel.grid=element_blank(),legend.position="bottom",
                                    plot.title=element_text(size=12))
  if (abs(sin(2*(a-b))) > 1e-12) {
    p <- p + geom_contour(data=projected,aes(x=s,y=t,z=relative_difference),
                          inherit.aes=FALSE,breaks=0,color="#333333",linewidth=0.3) +
      annotate("label",x=0,y=0,label="-",size=5,linewidth=0,
               fill="white",label.padding=unit(0.04,"lines")) +
      annotate("label",x=0,y=0.88,label="+",size=5,linewidth=0,
               fill="white",label.padding=unit(0.04,"lines"))
  }
  p
}
plots <- list(reference,make_projection(pi/8,pi/8,"Common origins: invisible"),
              make_projection(pi/4,0,"Crossed origins: visible"))
draw <- function() {
  grid.newpage();pushViewport(viewport(layout=grid.layout(1,3)))
  for(i in 1:3) print(plots[[i]],vp=viewport(layout.pos.row=1,layout.pos.col=i))
  popViewport()
}
pdf(file.path(out,"population_example.pdf"),width=9,height=3.7);draw();dev.off()
png(file.path(out,"population_example.png"),width=2700,height=1110,res=300);draw();dev.off()
write.csv(angular,file.path(out,"population_reference.csv"),row.names=FALSE)
cat("Population example rendered; no sampled or human results.\n")
