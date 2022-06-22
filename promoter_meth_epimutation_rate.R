######
getBinaryPromoterMethylation <- function(promoter_granges, cpg_granges){
  overlaps <- findOverlaps(cpg_granges, promoter_granges)
  cpgs_proms <- cpg_granges[unique(queryHits(overlaps))]
  cpgs_proms$tx_id <- promoter_granges$tx_id[subjectHits(overlaps)]
  
  cpgs_proms <- cpgs_proms[which(cpgs_proms$coverage >= 10)]
  
  human_cpgs_by_tx <- split(cpgs_proms, cpgs_proms$tx_id)
  meth_percentage <- sapply(human_cpgs_by_tx, function(xx) {
    xx <- xx[which(xx$meth_level <= 0.2 | xx$meth_level >= 0.8)]
    length(which(xx$meth_level >= 0.8))/length(xx$meth_level)
  })
  meth_percentage
}

getMeanPromoterMethylation <- function(promoter_granges, cpg_granges){
  overlaps <- findOverlaps(cpg_granges, promoter_granges)
  cpgs_proms <- cpg_granges[unique(queryHits(overlaps))]
  cpgs_proms$tx_id <- promoter_granges$tx_id[subjectHits(overlaps)]
  
  human_cpgs_by_tx <- split(cpgs_proms, cpgs_proms$tx_id)
  mean_meth <- sapply(human_cpgs_by_tx, function(xx) mean(xx$meth_level))
  mean_meth
}

proms_meth_human <- getBinaryPromoterMethylation(proms_good_power - 1500, cpgs_with_meth_level)
proms_meth_chimp <- getBinaryPromoterMethylation(proms_good_power - 1500, cpgs_with_meth_level_chimp)
proms_meth_rhesus <- getBinaryPromoterMethylation(proms_good_power - 1500, cpgs_with_meth_level_rhesus)

methylated_in_chimp <- proms_meth_chimp[which(proms_meth_chimp >= 0.8)]
methylated_in_chimp_and_rhesus <- proms_meth_rhesus[names(methylated_in_chimp)]
methylated_in_chimp_and_rhesus <- methylated_in_chimp_and_rhesus[which(methylated_in_chimp_and_rhesus >= 0.8)]
loeufs_to_use <- constraint[names(methylated_in_chimp_and_rhesus), "oe_lof_upper"]
methylated_in_chimp_and_rhesus <- methylated_in_chimp_and_rhesus[which(loeufs_to_use >= quantile(proms_good_power$oe_lof_upper, 0.75))]
length(which(proms_meth_human[names(methylated_in_chimp_and_rhesus)] <= 0.2))/length(proms_meth_human[names(methylated_in_chimp_and_rhesus)])

hypomethylated_in_chimp <- proms_meth_chimp[which(proms_meth_chimp <= 0.2)]
hypomethylated_in_chimp_and_rhesus <- proms_meth_rhesus[names(hypomethylated_in_chimp)]
hypomethylated_in_chimp_and_rhesus <- hypomethylated_in_chimp_and_rhesus[which(hypomethylated_in_chimp_and_rhesus <= 0.2)]
loeufs_to_use <- constraint[names(hypomethylated_in_chimp_and_rhesus), "oe_lof_upper"]
hypomethylated_in_chimp_and_rhesus <- hypomethylated_in_chimp_and_rhesus[which(loeufs_to_use >= quantile(proms_good_power$oe_lof_upper, 0.75))]
length(which(proms_meth_human[names(hypomethylated_in_chimp_and_rhesus)] >= 0.8))/length(proms_meth_human[names(hypomethylated_in_chimp_and_rhesus)])


quartz(file = "proximal_prom_comparative.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
hist(proms_meth_human[names(methylated_in_chimp_and_rhesus)], col = alpha("red", 0.6), lty = 0, 
     breaks = 35, freq = FALSE, xlab = "% methylated CpGs in human", cex.lab = 1,
     yaxt = 'n', xaxt = 'n', main = "proximal promoters methylated\nin chimp & rhesus", font.main = 1, cex.main = 1, xlim =  c(0, 1))
axis(1, at = c(0, 1), cex.axis = 1)
axis(2, at = c(0, 18), cex.axis = 1)

hist(proms_meth_human[names(hypomethylated_in_chimp_and_rhesus)], col = alpha("red", 0.6), lty = 0, 
     breaks = 35, freq = FALSE, xlab = "% methylated CpGs in human", cex.lab = 1,
     yaxt = 'n', xaxt = 'n', main = "proximal promoters hypomethylated\nin chimp & rhesus", font.main = 1, cex.main = 1, xlim =  c(0, 1))
axis(1, at = c(0, 1), cex.axis = 1)
axis(2, at = c(0, 18), cex.axis = 1)

dev.off()







####alternative analysis (not necessary - not used in the paper)
proms_meth_human <- getMeanPromoterMethylation(proms_good_power - 1500, cpgs_with_meth_level)
proms_meth_chimp <- getMeanPromoterMethylation(proms_good_power - 1500, cpgs_with_meth_level_chimp)
proms_meth_rhesus <- getMeanPromoterMethylation(proms_good_power - 1500, cpgs_with_meth_level_rhesus)

methylated_in_human <- names(proms_meth_human)[which(proms_meth_human >= 0.8)]
chimp_proms_methylated_in_human <- proms_meth_chimp[methylated_in_human]
rhesus_proms_methylated_in_human <- proms_meth_rhesus[names(chimp_proms_methylated_in_human)]
chimp_proms_methylated_in_human <- chimp_proms_methylated_in_human[names(rhesus_proms_methylated_in_human)]

length(which(rhesus_proms_methylated_in_human <= 0.5 & chimp_proms_methylated_in_human <= 0.5))

hypomethylated_in_human <- names(proms_meth_human_mean)[which(proms_meth_human_mean <= 0.2)]
chimp_proms_hypomethylated_in_human <- proms_meth_chimp_mean[hypomethylated_in_human]
rhesus_proms_hypomethylated_in_human <- proms_meth_rhesus_mean[names(chimp_proms_hypomethylated_in_human)]
chimp_proms_hypomethylated_in_human <- chimp_proms_hypomethylated_in_human[names(rhesus_proms_hypomethylated_in_human)]

length(which(rhesus_proms_hypomethylated_in_human >= 0.5 & chimp_proms_hypomethylated_in_human >= 0.5))


