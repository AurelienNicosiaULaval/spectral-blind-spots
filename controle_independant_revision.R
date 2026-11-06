library(stats)
# Check, by midpoint quadrature, that sin(2 theta) is orthogonal to
# every response factor used by the canonical views, including the
# sin/cos projections under each of the two origin choices.
theta <- 2*pi*(0:4095)/4096
res <- expand.grid(origin=c(0,pi/2),order=1:8)
res$integral <- vapply(seq_len(nrow(res)),function(i) mean(sin(2*theta)*cos(res$order[i]*(theta-res$origin[i]))),0.0)
stopifnot(max(abs(res$integral))<1e-14)
# Re-derive the original projected-cell list independently from the condition
# categories and verify the reported count and original local rejections.
s <- read.csv('resultats_simulation/rejection_summary.csv')
phase <- expand.grid(scenario=c('phase_2_moderate','phase_2_strong','phase_4','phase_6'),diagnostic=c('canonical','common_rotations'))
other <- data.frame(scenario=c('phase_4','phase_7_outside','product_2','local_rare','local_common'),diagnostic=c('crossed_pi4','crossed_band6','canonical','canonical','canonical'))
key <- rbind(phase,other); expected <- merge(key,s,by=c('scenario','diagnostic'))
stopifnot(nrow(expected)==26,all(expected$replications==200))
actual <- read.csv('resultats_simulation/population_blind_cells.csv')
actual <- actual[actual$diagnostic!='marginal_8',]
stopifnot(setequal(paste(actual$scenario,actual$diagnostic,actual$n),paste(expected$scenario,expected$diagnostic,expected$n)))
write.csv(res,'complement/controle_integrales_locales.csv',row.names=FALSE)
cat('Maximum absolute angular average:',max(abs(res$integral)),'\n26 projected scenario-size cells verified\n')
