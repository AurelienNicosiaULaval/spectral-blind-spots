# Paired validation on one reproducible simulated dataset.
# This illustration is separate from the 5,600-dataset study.
library(stats)
library(utils)
library(grDevices)
library(graphics)
script <- sub("^--file=", "", commandArgs()[startsWith(commandArgs(), "--file=")][1])
root <- dirname(normalizePath(script))
source(file.path(root,"diagnostics_core.R"))

# A uniform candidate loses an order-two dependence present in the reference law.
set.seed(2026100931)
definition <- list(kind="phase",order=2L,amplitude=0.8,is_null=FALSE)
x <- draw_x(720L,rho=0)
y <- draw_response(x,definition)
z <- runif(length(x),0,2*pi)
signs <- matrix(sample(c(-1,1),length(x)*199L,replace=TRUE),length(x),199L)
full <- sign_test(fourier_features(x,y,z),signs)
field <- field_features(x,y,z)
local <- sign_test(field$C,signs)
print(data.frame(diagnostic=c("Full Fourier","Periodic field"),
                 p_value=c(full$p,local$p),statistic=c(full$statistic,local$statistic)))

# One descriptive bandwidth, without significance labels for individual cells.
gx <- 2*pi*(0:15)/16
gt <- 2*pi*(0:31)/32
predictor <- wrapped_kernel(gx,x,h=0.35)
weights <- predictor/rowSums(predictor)
response_y <- wrapped_kernel(gt,y,h=0.15)
response_z <- wrapped_kernel(gt,z,h=0.15)
normalize_grid <- function(values) sweep(values,2,colSums(values)*2*pi/length(gt),"/")
difference <- weights%*%t(normalize_grid(response_z)-normalize_grid(response_y))
pdf(file.path(root,"figures","toy_field.pdf"),width=6,height=5)
closed_difference <- rbind(difference,difference[1,,drop=FALSE])
closed_difference <- cbind(closed_difference,closed_difference[,1])
limit <- max(abs(closed_difference))
filled.contour(c(gx,2*pi)*180/pi,c(gt,2*pi)*180/pi,closed_difference,
      xlab="Predictor (degrees)",ylab="Response (degrees)",
      main="Smoothed conditional difference",
      key.title=title(main="1/rad",cex.main=0.8),
      color.palette=colorRampPalette(c("#27649A","#F8F8F6","#B33D32")),
      levels=seq(-limit,limit,length.out=21L),
      plot.axes={axis(1,at=seq(0,360,90));axis(2,at=seq(0,360,90))})
dev.off()
# Supply the pairs so that the reader can change the diagnostic, not the reference law.
write.csv(data.frame(x,y,z),file.path(root,"figures","toy_pairs.csv"),row.names=FALSE)
