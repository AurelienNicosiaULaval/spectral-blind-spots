library(stats)
library(utils)
library(ggplot2)
source('diagnostics_core.R')
p<-read.csv('application/paires.csv');dates<-sort(unique(p$date))
train<-p[p$date%in%dates[1:50],];test<-p[p$date%in%dates[51:73],]
stopifnot(nrow(train)==600,nrow(test)==276,!any(train$date%in%test$date))
observed<-mean(exp(1i*(train$y-train$x)));mu<-Arg(observed);rbar<-Mod(observed)
ratio<-function(k)besselI(k,1,expon.scaled=TRUE)/besselI(k,0,expon.scaled=TRUE)
kappa<-uniroot(function(k)ratio(k)-rbar,c(1e-8,1000),tol=1e-12)$root
logi0<-function(k)log(besselI(k,0,expon.scaled=TRUE))+k
old<-read.csv('application/ajustement.csv');rad<-sqrt(old$k^2+(old$lambda*sin(test$x))^2)
ll_old<-old$k*cos(test$y-old$mu)+old$lambda*sin(test$x)*sin(test$y-old$mu)-log(2*pi)-logi0(rad)
ll_new<-kappa*cos(test$y-test$x-mu)-log(2*pi)-logi0(kappa)
write.csv(data.frame(kappa=kappa,shift=mu,resultant_train=rbar,n_train=600,n_test=276),'application/persistance_ajustement.csv',row.names=FALSE)
scores<-data.frame(date=test$date,log_original=ll_old,log_persistence=ll_new,difference=ll_new-ll_old)
write.csv(scores,'application/scores_paires.csv',row.names=FALSE)
write.csv(aggregate(scores[c('log_original','log_persistence','difference')],list(date=scores$date),mean),'application/scores_journaliers.csv',row.names=FALSE)
write.csv(data.frame(model=c('Original model','Persistence model'),mean_log_density=c(mean(ll_old),mean(ll_new)),n=276),'application/scores_synthese.csv',row.names=FALSE)
draw_vm<-function(x){
 out<-numeric(length(x));pending<-seq_along(x)
 while(length(pending)){
  proposal<-runif(length(pending),-pi,pi)
  keep<-runif(length(pending))<exp(kappa*(cos(proposal)-1))
  out[pending[keep]]<-wrap(x[pending[keep]]+mu+proposal[keep]);pending<-pending[!keep]
 }
 out
}
set.seed(2026110250);new<-replicate(50,draw_vm(test$x));old_draws<-as.matrix(read.csv('application/tirages.csv'))
stopifnot(identical(dim(new),dim(old_draws)),all(is.finite(new)))
write.csv(new,'application/tirages_persistance.csv',row.names=FALSE)
gx<-2*pi*(0:31)/32;gy<-2*pi*(0:63)/64
make_field<-function(z,hx,ht,model){
 K<-wrapped_kernel(gx,test$x,hx);w<-K/rowSums(K);neff<-1/rowSums(w^2)
 Ky<-wrapped_kernel(gy,test$y,ht);Ky<-sweep(Ky,2,colSums(Ky)*2*pi/length(gy),'/')
 f<-matrix(0,length(gx),length(gy))
 for(j in 1:50){Kz<-wrapped_kernel(gy,z[,j],ht);Kz<-sweep(Kz,2,colSums(Kz)*2*pi/length(gy),'/');f<-f+w%*%t(Kz-Ky)/50}
 d<-expand.grid(x=gx,y=gy);d$difference<-as.vector(f);d$n_eff<-rep(neff,length(gy));d$difference[d$n_eff<20]<-NA_real_;d$model<-model;d
}
main<-rbind(make_field(old_draws,.7,.3,'Original model'),make_field(new,.7,.3,'Persistence model'))
best<-read.csv('application/bandes_cv.csv');best<-best[best$selected,]
cv<-rbind(make_field(old_draws,best$hx,best$ht,'Original model'),make_field(new,best$hx,best$ht,'Persistence model'))
write.csv(main,'application/champs_comparaison.csv',row.names=FALSE);write.csv(cv,'application/champs_comparaison_cv.csv',row.names=FALSE)
plot_fields<-function(d,path){
 d$model<-factor(d$model,levels=c('Original model','Persistence model'))
 extrema<-do.call(rbind,lapply(split(d,d$model),function(v){r<-v[c(which.max(v$difference),which.min(v$difference)),];r$sign<-c('+','-');r}))
 lim<-max(abs(d$difference),na.rm=TRUE)
 dd<-rbind(d,transform(d[d$x==0,],x=2*pi));dd<-rbind(dd,transform(dd[dd$y==0,],y=2*pi))
 plot<-ggplot(dd,aes(x,y,fill=difference))+geom_raster()+
  geom_contour(aes(x=x,y=y,z=difference),inherit.aes=FALSE,breaks=0,color='black',linewidth=.3,na.rm=TRUE)+
  geom_label(data=extrema,aes(x=x,y=y,label=sign),inherit.aes=FALSE,fill='white',color='black',size=4,linewidth=0,label.padding=unit(.04,'lines'))+
  facet_wrap(~model,nrow=1)+scale_fill_gradient2(low='#2166ac',mid='white',high='#b2182b',limits=c(-lim,lim),na.value='grey80',name='Difference\nper radian')+
  coord_equal(xlim=c(0,2*pi),ylim=c(0,2*pi),expand=FALSE)+
  scale_x_continuous(breaks=c(0,pi,2*pi),labels=c('0',expression(pi),expression(2*pi)))+
  scale_y_continuous(breaks=c(0,pi,2*pi),labels=c('0',expression(pi),expression(2*pi)))+
  labs(x='Predictor direction (radians)',y='Response direction (radians)')+theme_minimal(base_size=11)
 ggsave(path,plot,width=9.2,height=4.2)
}
plot_fields(main,'application/champs_comparaison.pdf');plot_fields(cv,'application/champs_comparaison_cv.pdf')
# Empirical probabilities per fixed projected bin, candidate minus observed.
breaks<-seq(-1,1,length.out=25)
projected<-list();ix<-0L
for(model in c('Original model','Persistence model'))for(a in c(0,pi/2))for(b in c(0,pi/2)) {
 zz<-if(model=='Original model')old_draws else new
 bx<-findInterval(cos(test$x-a),breaks,all.inside=TRUE);by<-findInterval(cos(test$y-b),breaks,all.inside=TRUE)
 obs<-matrix(tabulate((by-1)*24+bx,nbins=24*24),24,24)/nrow(test)
 bz<-findInterval(as.vector(cos(zz-b)),breaks,all.inside=TRUE)
 candidate<-matrix(tabulate((bz-1)*24+rep(bx,50),nbins=24*24),24,24)/(50*nrow(test))
 d<-expand.grid(x=head(breaks,-1)+diff(breaks)/2,y=head(breaks,-1)+diff(breaks)/2)
 d$difference<-as.vector(candidate-obs);d$model<-model
 d$view<-paste0(if(a==0)'cos(X)' else 'sin(X)',', ',if(b==0)'cos(response)' else 'sin(response)')
 ix<-ix+1L;projected[[ix]]<-d
}
projected<-do.call(rbind,projected);write.csv(projected,'application/vues_projetees.csv',row.names=FALSE)
lim<-max(abs(projected$difference))
g<-ggplot(projected,aes(x,y,fill=difference))+geom_raster()+facet_grid(model~view)+
 scale_fill_gradient2(low='#2166ac',mid='white',high='#b2182b',limits=c(-lim,lim),name='Probability\ndifference')+
 coord_equal()+scale_x_continuous(breaks=c(-1,0,1))+scale_y_continuous(breaks=c(-1,0,1))+
 labs(x='Projected predictor',y='Projected response')+theme_minimal(base_size=10)
ggsave('application/vues_projetees.pdf',g,width=11,height=5.6)
cat('Wind model comparison, 25 CV pairs, descriptive scores and projected displays generated.\n')
