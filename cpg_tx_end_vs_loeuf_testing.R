###
edb <- EnsDb.Hsapiens.v75
tx_all <- transcripts(edb, columns = c("tx_id", "gene_id"))

#tx_all <- tx_all[proms_good_power$tx_id] #only for transcripts which have reliable promoter annotation
tx_all <- tx_all[names(exons_by_tx)] #this is for all canonical transcripts
tx_all <- tx_all+500
tx_end <- resize(tx_all, 1000, fix = "end")
seqlevelsStyle(tx_end) <- "ucsc"
genome(tx_end) <- "hg19"
seqlevels(tx_end) <- seqlevels(tx_end)[1:22]


cpgs_to_use <- cpgs_with_meth_level[which(cpgs_with_meth_level$coverage >= 10)]
cpgs_to_use <- cpgs_to_use[-queryHits(findOverlaps(cpgs_to_use, proms_all))] 
overlaps <- findOverlaps(tx_end, cpgs_to_use)
cpgs_by_tx_end <- split(cpgs_to_use[subjectHits(overlaps)], queryHits(overlaps))
names(cpgs_by_tx_end) <- names(tx_end)[unique(queryHits(overlaps))]

human_loeufs_2 <- unlist(sapply(names(cpgs_by_tx_end), function(xx) constraint$oe_lof_upper[which(constraint$transcript == xx)]))
indices <- which(lengths(cpgs_by_tx_end) >= 10)
cpgs_by_tx_end <- cpgs_by_tx_end[indices]
human_loeufs_2 <- human_loeufs_2[indices]
meth_percentage_tx_end <- sapply(cpgs_by_tx_end, function(xx) {
  xx <- xx[which(xx$meth_level <= 0.2 | xx$meth_level >= 0.8)]
  length(which(xx$meth_level >= 0.8))/length(xx$meth_level)
})
meth_percentage_tx_end <- meth_percentage_tx_end[order(human_loeufs_2)]
exp_lof_2 <- unlist(sapply(names(meth_percentage_tx_end), function(xx) constraint$exp_lof[which(constraint$transcript == xx)]))
meth_percentage_tx_end <- meth_percentage_tx_end[which(exp_lof_2 >= 10)]

mean_meth_tx_end <- sapply(cpgs_by_tx_end, function(xx) mean(xx$meth_level))
mean_meth_tx_end <- mean_meth_tx_end[order(human_loeufs_2)]
mean_meth_tx_end <- mean_meth_tx_end[which(exp_lof_2 >= 10)]


meth_vs_loeuf_vec_tx_end <- c(getPercentOfMethPromoters(meth_percentage_tx_end, human_loeufs_2, 0.1, "less", 0.8), 
                       getPercentOfMethPromoters(meth_percentage_tx_end, human_loeufs_2, c(0.1, 0.2), "in between", 0.8),
                       getPercentOfMethPromoters(meth_percentage_tx_end, human_loeufs_2, c(0.2, 0.3), "in between", 0.8),
                       getPercentOfMethPromoters(meth_percentage_tx_end, human_loeufs_2, c(0.3, 0.4), "in between", 0.8),
                       getPercentOfMethPromoters(meth_percentage_tx_end, human_loeufs_2, c(0.4, 0.5), "in between", 0.8),
                       getPercentOfMethPromoters(meth_percentage_tx_end, human_loeufs_2, c(0.5, 0.6), "in between", 0.8),
                       getPercentOfMethPromoters(meth_percentage_tx_end, human_loeufs_2, c(0.6, 0.7), "in between", 0.8),
                       getPercentOfMethPromoters(meth_percentage_tx_end, human_loeufs_2, c(0.7, 0.8), "in between", 0.8),
                       getPercentOfMethPromoters(meth_percentage_tx_end, human_loeufs_2, c(0.8, 0.9), "in between", 0.8),
                       getPercentOfMethPromoters(meth_percentage_tx_end, human_loeufs_2, 0.9, "greater", 0.8))


###
quartz(file = "methylation_vs_loeuf_tx_end.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(meth_vs_loeuf_vec_tx_end, pch = 19, cex = 1.2, col = alpha("red", 0.75), bty = 'l', ylim = c(0.6, 0.8),
     main = "transcriptional end region", font.main = 1, ylab = "% genes w/ methylation", xlab = "LOEUF decile", 
     xaxt = 'n', yaxt = 'n')
axis(1, at = 1:10, cex.axis = 0.7)
axis(2, at = c(0.6, 0.8))
dev.off()
###




###
quartz(file = "tx_end_methylation_histogram.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
hist(meth_percentage_tx_end, col = alpha("red", 0.6), lty = 0, 
     breaks = 35, freq = FALSE, xlab = "% methylated CpGs", cex.lab = 1.25,
     yaxt = 'n', xaxt = 'n', main = "transcriptional end site\n(+/- 500bp)", font.main = 1, cex.main = 1.25, xlim =  c(0, 1))
axis(1, at = c(0, 1), cex.axis = 1.2)
axis(2, at = c(0, 15), cex.axis = 1.2)
#abline(v = c(40, 80), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()

