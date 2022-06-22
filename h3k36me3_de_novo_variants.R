###relationship with h3k36me3
exons_h3k36me3_overlapping_interval <- intersect(exons_to_use, h3k36me3, ignore.strand = TRUE) + 250
#next line ensures that the addition of 250bp on either side does not include introns
exons_h3k36me3_overlapping_interval <- intersect(exons_h3k36me3_overlapping_interval, exons_to_use, ignore.strand = TRUE) 
exons_non_h3k36me3_overlapping_interval <- setdiff(exons_to_use, exons_h3k36me3_overlapping_interval, ignore.strand = TRUE)

h3m36me3_mut_rate <- length(unique(queryHits(findOverlaps(de_novo_mutations, 
                                                          exons_h3k36me3_overlapping_interval))))/sum(width(exons_h3k36me3_overlapping_interval))
non_h3m36me3_mut_rate <- length(unique(queryHits(findOverlaps(de_novo_mutations, 
                                                              exons_non_h3k36me3_overlapping_interval))))/sum(width(exons_non_h3k36me3_overlapping_interval))
###
permutations <- replicate(10000, {
  indices_all <- sample(1:length(exons_to_use), 
                        c(length(exons_h3k36me3_overlapping_interval), length(exons_non_h3k36me3_overlapping_interval)))
  indices_1 <- indices_all[1:length(exons_h3k36me3_overlapping_interval)]
  indices_2 <- indices_all[-indices_1]
  permuted_mut_rate_1 <- length(unique(queryHits(findOverlaps(de_novo_mutations_in_tx, 
                                                              exons_to_use[indices_1]))))/sum(width(exons_to_use[indices_1]))
  permuted_mut_rate_2 <- length(unique(queryHits(findOverlaps(de_novo_mutations_in_tx, 
                                                              exons_to_use[indices_2]))))/sum(width(exons_to_use[indices_2]))
  permuted_mut_rate_2 - permuted_mut_rate_1
})


quartz(file = "h3k36me3_mut_rate.pdf", width = 4.4, height = 2.2, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
hist(permutations, breaks = 40, xlim = c(min(permutations), max(non_h3m36me3_mut_rate - h3m36me3_mut_rate, max(permutations))), 
     main = "germline H3K36me3", xlab = "mut rate difference\n(peaks vs outside peaks)", font.main = 1)
abline(v = non_h3m36me3_mut_rate - h3m36me3_mut_rate, col = alpha("red", 0.57), lwd = 2.5)

hist(permutations_h1esc, breaks = 40, 
     main = "ESC H3K36me3", xlab = "mut rate difference\n(peaks vs outside peaks)", xlim = c(-2e-05, 7e-05), font.main = 1)
abline(v = non_h3m36me3_mut_rate_h1esc - h3m36me3_mut_rate_h1esc, col = alpha("red", 0.57), lwd = 2.5)
dev.off()


###de novo mut rate vs loeuf excluding h3k36me3 regions
#
de_novo_overlaps2 <- findOverlaps(exons_non_h3k36me3_overlapping_interval, de_novo_mutations)
exons_non_h3k36me3_overlapping_interval$de_novo_count <- 0
mut_by_exon2 <- split(de_novo_mutations[subjectHits(de_novo_overlaps2)], queryHits(de_novo_overlaps2))
exons_non_h3k36me3_overlapping_interval$de_novo_count[unique(queryHits(de_novo_overlaps2))] <- lengths(mut_by_exon2)

tx_all <- transcripts(edb, columns = c("tx_id", "gene_id"))
tx_all <- tx_all[names(exons_by_tx)]
seqlevelsStyle(tx_all) <- "ucsc"
genome(tx_all) <- "hg19"
seqlevels(tx_end) <- seqlevels(tx_end)[1:22]
overlaps <- findOverlaps(exons_non_h3k36me3_overlapping_interval, tx_all)
duplicated_indices <- queryHits(overlaps)[which(duplicated(queryHits(overlaps)))]
exons_non_h3k36me3_overlapping_interval <- exons_non_h3k36me3_overlapping_interval[-duplicated_indices]

overlaps <- findOverlaps(exons_non_h3k36me3_overlapping_interval, tx_all)
exons_to_use_2 <- exons_non_h3k36me3_overlapping_interval[unique(queryHits(overlaps))]
exons_to_use_2$tx_id <-  tx_all$tx_id[subjectHits(overlaps)]
names(exons_to_use_2) <- exons_to_use_2$tx_id

exons_by_tx_2 <- split(exons_to_use_2, exons_to_use_2$tx_id)
coding_length_2 <- sapply(exons_by_tx_2, function(xx) sum(width(xx)))
oe_lof_upper_2 <- constraint[names(exons_by_tx_2), "oe_lof_upper"]
de_novo_count_2 <- sapply(exons_by_tx_2, function(xx) sum(xx$de_novo_count))

de_novo_count_df_2 <- data.frame(de_novo_count = de_novo_count_2, cds_length = coding_length_2, 
                               oe_lof_upper = oe_lof_upper_2, tx_id = names(exons_by_tx_2),
                               stringsAsFactors = FALSE)
de_novo_count_df_2$mut_rate <- NA
de_novo_count_df_2$mut_rate[which(de_novo_count_df_2$de_novo_count > 0)] <- de_novo_count_df_2$de_novo_count[
  which(de_novo_count_df_2$de_novo_count > 0)]/de_novo_count_df_2$cds_length[which(de_novo_count_df_2$de_novo_count > 0)]

de_novo_count_df_2$mut_rate_2 <- (de_novo_count_df_2$de_novo_count + 1)/de_novo_count_df_2$cds_length

#de_novo_count_df_2 <- de_novo_count_df_2[-which(is.na(de_novo_count_df_2$oe_lof_upper)), ]

##restrict to genes with k36me3 peak at coding region
overlaps <- findOverlaps(exons_to_use, h3k36me3)
tx_ids <- unique(exons_to_use$tx_id[queryHits(overlaps)])

df1 <- de_novo_count_df[which(de_novo_count_df$tx_id %in% tx_ids), ]
df2 <- de_novo_count_df_2[which(de_novo_count_df_2$tx_id %in% tx_ids), ]

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

#make the same plot this time using the mutation rate estimate with the pseudocount
diff2 <- unlist(sapply(df2$tx_id, function(xx) df1$mut_rate_2[which(df1$tx_id == xx)])) - df2$mut_rate_2

quartz(file = "coding_mut_rate_difference_h3k36me3_2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
hist(diff2, col = alpha("red", 0.6), lty = 0, 
     breaks = 2000, freq = FALSE, xlab = "coding mut rate difference\n(with minus w/out H3K36me3)", cex.lab = 1,
     yaxt = 'n', xaxt = 'n', main = "", font.main = 1, xlim =  c(-0.005, 0.005))
axis(1, at = c(-0.003, 0, 0.003), labels = c("-0.003", "0", "0.003"), cex.axis = 0.8)
axis(2, at = c(0, 2200), cex.axis = 0.8)
#abline(v = 0, lty = "longdash", col = rgb(0,0,0,0.7), lwd = 2.5)
dev.off()





overlaps <- findOverlaps(exons_non_h3k36me3_overlapping_interval, exons_to_use)
exons_to_use_2 <- exons_non_h3k36me3_overlapping_interval[queryHits(overlaps)]
exons_to_use_2$tx_id <-  exons_to_use$tx_id[subjectHits(overlaps)]
names(exons_to_use_2) <- exons_to_use_2$tx_id

exons_by_tx_2 <- split(exons_to_use_2, exons_to_use_2$tx_id)
coding_length_2 <- sapply(exons_by_tx_2, function(xx) sum(width(xx)))
oe_lof_upper_2 <- constraint[names(exons_by_tx_2), "oe_lof_upper"]

de_novo_overlaps <- findOverlaps(de_novo_mutations, exons_to_use_2)
de_novo_mutations_in_tx <- de_novo_mutations[queryHits(de_novo_overlaps)]
de_novo_mutations_in_tx$tx_id <- exons_to_use_2$tx_id[subjectHits(de_novo_overlaps)]
de_novo_mutations_in_tx <- unique(de_novo_mutations_in_tx)

de_novo_mutations_by_tx <- split(de_novo_mutations_in_tx, de_novo_mutations_in_tx$tx_id)

indices_with_mut <- which(names(exons_by_tx_2) %in% names(de_novo_mutations_by_tx))
de_novo_mutations_by_tx <- de_novo_mutations_by_tx[names(exons_by_tx_2)[indices_with_mut]]

de_novo_count_2 <- rep(0, length(names(exons_by_tx_2)))
de_novo_count_2[indices_with_mut] <- lengths(de_novo_mutations_by_tx)

de_novo_count_df_2 <- data.frame(de_novo_count = de_novo_count_2, cds_length = coding_length_2, 
                               oe_lof_upper = oe_lof_upper_2, tx_id = names(exons_by_tx_2),
                               stringsAsFactors = FALSE)
de_novo_count_df_2$mut_rate <- NA
de_novo_count_df_2$mut_rate[which(de_novo_count_df_2$de_novo_count > 0)] <- de_novo_count_df_2$de_novo_count[
  which(de_novo_count_df_2$de_novo_count > 0)]/de_novo_count_df_2$cds_length[which(de_novo_count_df_2$de_novo_count > 0)]

de_novo_count_df_2$mut_rate_2 <- (de_novo_count_df_2$de_novo_count + 0.5)/de_novo_count_df_2$cds_length

de_novo_count_df_2 <- de_novo_count_df_2[-which(is.na(de_novo_count_df_2$oe_lof_upper)), ]

mut_rate_without_k36me3 <- de_novo_count_df_2$mut_rate_2
mut_rate_with_k36me3 <- sapply(de_novo_count_df_2$tx_id, function(xx) de_novo_count_df$mut_rate_2[which(de_novo_count_df$tx_id == xx)])

diff <- mut_rate_without_k36me3 - mut_rate_with_k36me3




