library(stats)
library(utils)
library(grDevices)
library(ggplot2)
script <- sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")][1])
root <- dirname(normalizePath(script))
results_dir <- file.path(root,"resultats_simulation")
figure_dir <- file.path(root,"figures")
dir.create(figure_dir,showWarnings=FALSE)
completion <- readRDS(file.path(results_dir,"completion.rds"))
raw <- read.csv(file.path(results_dir,"all_results.csv"))
stopifnot(completion$datasets==5600L,nrow(raw)==56000L,
          !anyDuplicated(raw[c("condition","replication","diagnostic")]))

wilson <- function(successes,n) {
  z <- qnorm(0.975); p <- successes/n
  center <- (p+z^2/(2*n))/(1+z^2/n)
  width <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n)
  c(lower=center-width,upper=center+width)
}
groups <- split(raw,interaction(raw$condition,raw$diagnostic,drop=TRUE))
summary <- do.call(rbind,lapply(groups,function(g) {
  R <- nrow(g); rejected <- sum(g$reject); ci <- wilson(rejected,R)
  data.frame(condition=g$condition[1],scenario=g$scenario[1],n=g$n[1],
    diagnostic=g$diagnostic[1],replications=R,rejections=rejected,
    rate=rejected/R,lower=ci[1],upper=ci[2],
    mc_standard_error=sqrt((rejected/R)*(1-rejected/R)/R),
    median_seconds=median(g$elapsed_seconds),
    median_components=median(g$components),
    min_field_positions=min(g$field_min_positions),
    max_field_positions=max(g$field_max_positions))
}))
summary <- summary[order(summary$condition,summary$diagnostic),]
write.csv(summary,file.path(results_dir,"rejection_summary.csv"),row.names=FALSE)

differences <- do.call(rbind,lapply(split(raw,raw$condition),function(g) {
  baseline <- g[g$diagnostic=="canonical",c("replication","reject")]
  do.call(rbind,lapply(setdiff(unique(g$diagnostic),"canonical"),function(name) {
    target <- g[g$diagnostic==name,c("replication","reject")]
    aligned <- merge(baseline,target,by="replication",suffixes=c("_baseline","_target"))
    d <- as.numeric(aligned$reject_target)-as.numeric(aligned$reject_baseline)
    data.frame(condition=g$condition[1],scenario=g$scenario[1],n=g$n[1],
      diagnostic=name,baseline="canonical",replications=length(d),
      difference=mean(d),mc_standard_error=sd(d)/sqrt(length(d)),
      only_target=sum(d==1),only_baseline=sum(d==-1))
  }))
}))
write.csv(differences,file.path(results_dir,"paired_differences.csv"),row.names=FALSE)

labels <- c(canonical="Canonical views",common_rotations="Common origins",
  crossed_pi4="Crossed origins: pi/4",crossed_band6="Crossed origins: pi/7",
  crossed_irrational="Crossed irrational",fourier_8="Full Fourier",
  marginal_8="Marginal Fourier",periodic_field="Periodic field",
  paired_mmd="Quadratic MMD",coverage_known_model="Known-law coverage")
order <- names(labels)
scenarios <- read.csv(file.path(results_dir,"conditions.csv"))
source(file.path(root,"diagnostics_core.R"))
laws <- scenarios[!duplicated(scenarios$scenario),]
population_metrics <- do.call(rbind,lapply(seq_len(nrow(laws)),function(index) {
  law <- laws[index,]
  if(law$is_null) {
    tv <- 0; max_conditional <- 0; uniform_average <- 0
  } else if(law$kind=="phase") {
    tv <- law$amplitude/pi; max_conditional <- tv; uniform_average <- tv
  } else if(law$kind=="product") {
    tv <- 2*law$amplitude/pi^2; max_conditional <- law$amplitude/pi; uniform_average <- tv
  } else {
    center <- if(law$kind=="local_rare") pi else 0
    tv <- integrate(function(x) (1+law$rho*cos(x))/(2*pi)*law$amplitude*bump(x,center)/pi,
                    0,2*pi,subdivisions=1000,abs.tol=1e-10,rel.tol=1e-10)$value
    max_conditional <- law$amplitude/pi
    uniform_average <- law$amplitude*(16*0.9/15)/(2*pi^2)
  }
  data.frame(scenario=law$scenario,joint_tv=tv,max_conditional_tv=max_conditional,
    angular_average_conditional_tv=uniform_average,
    continuum_min_band_attenuation=if(law$kind=="phase")
      exp(-law$order^2*(0.35^2+0.15^2)/2) else NA_real_,
    # Weight of the (k,-k) Fourier mode in the squared population MMD of the product kernel with s=0.5,
    # relative to order one: {I_k(4)/I_1(4)}^2 (revision of 10 October 2026).
    mmd_relative_spectral_weight=if(law$kind=="phase")
      (besselI(4,law$order,expon.scaled=TRUE)/besselI(4,1,expon.scaled=TRUE))^2 else NA_real_)
}))
write.csv(population_metrics,file.path(results_dir,"population_metrics.csv"),row.names=FALSE)

# Cells with zero population contrast for the projected or marginal features (revision of 10 October 2026).
# These are alternatives, not complete-fidelity nulls; rates are reported, not interpreted as exact levels.
blind <- rbind(
  expand.grid(scenario=c("phase_2_moderate","phase_2_strong","phase_4","phase_6"),
              diagnostic=c("canonical","common_rotations"),stringsAsFactors=FALSE),
  data.frame(scenario=c("phase_4","phase_7_outside","product_2","local_rare","local_common"),
             diagnostic=c("crossed_pi4","crossed_band6","canonical","canonical","canonical")),
  expand.grid(scenario=c("odd_1","phase_2_moderate","phase_2_strong","phase_4","phase_6",
                         "phase_7_outside","product_2"),diagnostic="marginal_8",stringsAsFactors=FALSE))
blind <- merge(blind,summary[c("scenario","diagnostic","n","replications","rejections","rate","lower","upper")],
               by=c("scenario","diagnostic"))
blind <- blind[order(blind$diagnostic,blind$scenario,blind$n),]
write.csv(blind,file.path(results_dir,"population_blind_cells.csv"),row.names=FALSE)
scenario_order <- unique(scenarios$scenario)
scenario_labels <- c(null_uniform="Uniform fidelity",null_conditional="Conditional fidelity",
  odd_1="Phase 1, amplitude 0.4",phase_2_moderate="Phase 2, amplitude 0.4",
  phase_2_strong="Phase 2, amplitude 0.8",phase_4="Phase 4, amplitude 0.8",
  phase_6="Phase 6, amplitude 0.8",phase_7_outside="Phase 7, amplitude 0.8",
  product_2="Even-sine product",local_rare="Local rare zone",
  local_common="Local common zone")
summary$label <- factor(labels[summary$diagnostic],levels=unname(labels))
summary$scenario_label <- factor(scenario_labels[summary$scenario],
                                 levels=rev(unname(scenario_labels[scenario_order])))
heat <- ggplot(summary,aes(label,scenario_label,fill=rate))+
  geom_tile(color="white",linewidth=0.4)+
  geom_text(aes(label=sprintf("%.2f",rate),color=rate>0.6),size=3.8)+
  scale_color_manual(values=c("FALSE"="#1F2933","TRUE"="white"),guide="none")+
  scale_fill_gradient(low="#F6F7F5",high="#265F8A",limits=c(0,1),name="Rejection rate")+
  facet_wrap(~n,ncol=1,labeller=label_both)+
  labs(x=NULL,y=NULL)+theme_minimal(base_size=11)+
  theme(panel.grid=element_blank(),axis.text.x=element_text(angle=60,hjust=1),
        legend.position="bottom",strip.text=element_text(size=12))
ggsave(file.path(figure_dir,"rejection_rates.pdf"),heat,width=8.2,height=8.3)
ggsave(file.path(figure_dir,"rejection_rates.png"),heat,width=8.2,height=8.3,dpi=180)

null <- summary[grepl("^null",summary$scenario),]
calibration <- ggplot(null,aes(rate,label,color=factor(n),shape=factor(n)))+
  geom_vline(xintercept=0.05,linetype=2,color="#555555")+
  geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=0.4,
                position=position_dodge(width=0.55))+
  geom_point(position=position_dodge(width=0.55),size=2)+
  facet_wrap(~scenario,ncol=2,labeller=as_labeller(c(null_uniform="Uniform fidelity",null_conditional="Conditional fidelity")))+
  scale_color_manual(values=c("#27649A","#AC3A32"),name="Sample size")+
  scale_shape_manual(values=c(16,17),name="Sample size")+
  labs(x="Rejection rate and 95% Wilson interval",y=NULL)+theme_minimal(base_size=11)+
  theme(panel.grid.major.y=element_blank(),legend.position="bottom")
ggsave(file.path(figure_dir,"calibration.pdf"),calibration,width=8,height=5)
ggsave(file.path(figure_dir,"calibration.png"),calibration,width=8,height=5,dpi=180)

conditioning <- read.csv(file.path(root,"resultats_geometry","conditioning.csv"))
conditioning <- conditioning[conditioning$frequency<=16,]
conditioning$minimum_singular_value[conditioning$minimum_singular_value<1e-10] <- 0
# One panel per angle (revision of 10 October 2026): overlaid lines were hard to separate.
angle_labels <- c(canonical="Canonical, pi/2",quarter_turn="Crossed, pi/4",
                  band6="Crossed, pi/7",irrational="Crossed, pi/sqrt(2)")
conditioning$angle_label <- factor(angle_labels[conditioning$angle],levels=unname(angle_labels))
conditioning$resonant <- conditioning$minimum_singular_value==0
plot_condition <- ggplot(conditioning,aes(frequency,minimum_singular_value))+
  geom_hline(yintercept=0,color="#BBBBBB",linewidth=0.3)+
  geom_segment(aes(xend=frequency,yend=0),color="#7A8B99",linewidth=0.4)+
  geom_point(aes(shape=resonant,color=resonant),size=2)+
  scale_shape_manual(values=c("FALSE"=16,"TRUE"=4),labels=c("FALSE"="Nonzero","TRUE"="Exact resonance"),name=NULL)+
  scale_color_manual(values=c("FALSE"="#27649A","TRUE"="#AC3A32"),labels=c("FALSE"="Nonzero","TRUE"="Exact resonance"),name=NULL)+
  scale_x_continuous(breaks=c(1,4,7,10,13,16))+
  facet_wrap(~angle_label,nrow=2)+
  labs(x="Fourier order k",y="Smallest singular value")+
  theme_minimal(base_size=12)+theme(legend.position="bottom",panel.grid.minor=element_blank())
ggsave(file.path(figure_dir,"conditioning.pdf"),plot_condition,width=7.5,height=4.2)
ggsave(file.path(figure_dir,"conditioning.png"),plot_condition,width=7.5,height=4.2,dpi=180)

# Table courte, résultats des mêmes fichiers que les figures.
selected <- summary[summary$n==720 & summary$diagnostic %in%
                     c("canonical","common_rotations","crossed_band6","fourier_8",
                       "periodic_field","paired_mmd"),]
wide <- reshape(selected[c("scenario","diagnostic","rate")],idvar="scenario",
                timevar="diagnostic",direction="wide")
wide <- wide[match(scenario_order,wide$scenario),]
names(wide) <- sub("rate\\.","",names(wide))
write.csv(wide,file.path(results_dir,"main_table.csv"),row.names=FALSE)
# The recorded R_citation.txt and ggplot2_citation.txt document the original execution (R 4.5.0,
# ggplot2 4.0.2) and are kept; citations of the current rebuild are written separately (revision of 10 October 2026).
writeLines(capture.output(citation()),file.path(results_dir,"R_citation_rebuild.txt"))
writeLines(capture.output(citation("ggplot2")),file.path(results_dir,"ggplot2_citation_rebuild.txt"))
cat("Report generated from all 5600 datasets; no additional simulations run.\n")
