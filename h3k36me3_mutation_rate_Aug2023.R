exons_h3k36me3_overlapping_interval <- intersect(exons_to_use, h3k36me3, ignore.strand = TRUE) + 250 #adding 250bp because recruitment of repair machinery can also affect flanking sites
#next line ensures that the addition of 250bp on either side does not include introns
exons_h3k36me3_overlapping_interval <- intersect(exons_h3k36me3_overlapping_interval, exons_to_use, ignore.strand = TRUE) 
exons_non_h3k36me3_overlapping_interval <- setdiff(exons_to_use, exons_h3k36me3_overlapping_interval, ignore.strand = TRUE)

h3m36me3_mut_rate <- length(unique(queryHits(findOverlaps(de_novo_mutations, 
                                                          exons_h3k36me3_overlapping_interval))))/sum(width(exons_h3k36me3_overlapping_interval))
non_h3m36me3_mut_rate <- length(unique(queryHits(findOverlaps(de_novo_mutations, 
                                                              exons_non_h3k36me3_overlapping_interval))))/sum(width(exons_non_h3k36me3_overlapping_interval))

genome(exons_non_h3k36me3_overlapping_interval) <- "hg19" #based on paper it is hg19
seqlevelsStyle(exons_non_h3k36me3_overlapping_interval) <- "ucsc"
overlaps <- findOverlaps(exons_non_h3k36me3_overlapping_interval, tx_all)
duplicated_indices <- queryHits(overlaps)[which(duplicated(queryHits(overlaps)))]
exons_non_h3k36me3_overlapping_interval <- exons_non_h3k36me3_overlapping_interval[-duplicated_indices]

overlaps <- findOverlaps(exons_non_h3k36me3_overlapping_interval, tx_all)
exons_to_use_2 <- exons_non_h3k36me3_overlapping_interval[unique(queryHits(overlaps))]
exons_to_use_2$tx_id <-  tx_all$tx_id[subjectHits(overlaps)]
names(exons_to_use_2) <- exons_to_use_2$tx_id

mut_rate_df_syn_no_introns_non_h3k36me3 <- getMutationRateDF(de_novo_df, "synonymous SNV", exons_to_use_2, 
                                                             region_type = "exon", use_pseudocount = FALSE)


##restrict to genes with k36me3 peak at coding region
overlaps <- findOverlaps(exons_to_use, h3k36me3)
ids_with_h3k36me3 <- unique(exons_to_use$tx_id[queryHits(overlaps)])

df1 <- mut_rate_df_syn_no_introns[which(mut_rate_df_syn_no_introns$tx_id %in% ids_with_h3k36me3), ]
df2 <- mut_rate_df_syn_no_introns_non_h3k36me3[which(mut_rate_df_syn_no_introns_non_h3k36me3$tx_id %in% 
                                                       ids_with_h3k36me3), ]

mut_rate_without_k36me3 <- df2$mut_rate
mut_rate_with_k36me3 <- unlist(sapply(df2$tx_id, function(xx) df1$mut_rate[which(df1$tx_id == xx)]))

diff <- mut_rate_with_k36me3 - mut_rate_without_k36me3

quartz(file = "coding_mut_rate_difference_h3k36me3.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
hist(diff, col = alpha("red", 0.6), lty = 0, 
     breaks = 400, freq = FALSE, xlab = "coding mut rate difference\n(with minus w/out H3K36me3)", cex.lab = 1,
     yaxt = 'n', xaxt = 'n', main = "", font.main = 1, xlim =  c(-0.0015, 0.0015))
axis(1, at = c(-0.0015, 0, 0.0015), labels = c("-0.0015", "0", "0.0015"), cex.axis = 0.8)
axis(2, at = c(0, 5700), cex.axis = 0.8)
#abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7), lwd = 2.5)
dev.off()



###now make the same plot but this time calculating the mutation rate with pseudocount
mut_rate_df_syn_no_introns_non_h3k36me3 <- getMutationRateDF(de_novo_df, "synonymous SNV", exons_to_use_2, 
                                                             region_type = "exon", use_pseudocount = TRUE, 
                                                             pseudocount = 1)

overlaps <- findOverlaps(exons_to_use, h3k36me3)
ids_with_h3k36me3 <- unique(exons_to_use$tx_id[queryHits(overlaps)])

df1 <- mut_rate_df_syn_no_introns[which(mut_rate_df_syn_no_introns$tx_id %in% ids_with_h3k36me3), ]
df2 <- mut_rate_df_syn_no_introns_non_h3k36me3[which(mut_rate_df_syn_no_introns_non_h3k36me3$tx_id %in% 
                                                       ids_with_h3k36me3), ]

mut_rate_without_k36me3 <- df2$mut_rate_2
mut_rate_with_k36me3 <- unlist(sapply(df2$tx_id, function(xx) df1$mut_rate_2[which(df1$tx_id == xx)]))

diff2 <- mut_rate_with_k36me3 - mut_rate_without_k36me3

quartz(file = "coding_mut_rate_difference_h3k36me3_2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
hist(diff2, col = alpha("red", 0.6), lty = 0, 
     breaks = 4000, freq = FALSE, xlab = "coding mut rate difference\n(with minus w/out H3K36me3)", cex.lab = 1,
     yaxt = 'n', xaxt = 'n', main = "", font.main = 1, xlim =  c(-0.005, 0.005))
axis(1, at = c(-0.003, 0, 0.003), labels = c("-0.003", "0", "0.003"), cex.axis = 0.8)
axis(2, at = c(0, 3000), cex.axis = 0.8)
#abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7), lwd = 2.5)
dev.off()




