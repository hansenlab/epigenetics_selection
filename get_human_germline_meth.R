library(rtracklayer)


#####
cpgs_with_meth_level <- import("~/Downloads/germline_methylation/human_germline_meth.bw", 
                               format = "BigWig")
cpgs_with_coverage <- import("~/Downloads/germline_methylation/human_germline_meth_coverage.bw", 
                             format = "BigWig")

cpgs_with_meth_level$meth_level <- cpgs_with_meth_level$score
cpgs_with_meth_level$coverage <- cpgs_with_coverage$score

load(file = "proms_good_power_exp_lof_more_than_10.rda")
overlaps <- findOverlaps(cpgs_with_meth_level, proms_good_power-1500)
cpgs_proms_good_power <- cpgs_with_meth_level[unique(queryHits(overlaps))]
cpgs_proms_good_power$tx_id <- proms_good_power$tx_id[subjectHits(overlaps)]

#the following line gives only the cpgs at the promoter boundary
#cpgs_proms_good_power <- cpgs_proms_good_power[queryHits(findOverlaps(cpgs_proms_good_power, setdiff(proms_good_power, proms_good_power - 1500)))]
#

cpgs_proms_good_power <- cpgs_proms_good_power[which(cpgs_proms_good_power$coverage >= 10)]

human_cpgs_by_tx <- split(cpgs_proms_good_power, cpgs_proms_good_power$tx_id)
human_loeufs <- sapply(names(human_cpgs_by_tx), function(xx) proms_good_power$oe_lof_upper[which(proms_good_power$tx_id == xx)])
human_loeufs <- human_loeufs[order(human_loeufs)]
human_cpgs_by_tx <- human_cpgs_by_tx[names(human_loeufs)]
indices <- which(lengths(human_cpgs_by_tx) >= 10)
human_cpgs_by_tx <- human_cpgs_by_tx[indices]
human_loeufs <- human_loeufs[indices]
mean_meth <- sapply(human_cpgs_by_tx, function(xx) mean(xx$meth_level))
meth_percentage <- sapply(human_cpgs_by_tx, function(xx) {
  xx <- xx[which(xx$meth_level <= 0.2 | xx$meth_level >= 0.8)]
  length(which(xx$meth_level >= 0.8))/length(xx$meth_level)
})

meth_percentage <- meth_percentage[order(human_loeufs)]
mean_meth <- mean_meth[order(human_loeufs)]



######
quartz(file = "core_promoter_methylation_histogram_human.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
hist(meth_percentage, col = alpha("red", 0.6), lty = 0, 
     breaks = 35, freq = FALSE, xlab = "% methylated CpGs", cex.lab = 1.25,
     yaxt = 'n', xaxt = 'n', main = "proximal promoter\n(+/- 500bp from TSS)", font.main = 1, cex.main = 1.25, xlim =  c(0, 1))
axis(1, at = c(0, 1), cex.axis = 1.2)
axis(2, at = c(0, 15), cex.axis = 1.2)
#abline(v = c(40, 80), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()













