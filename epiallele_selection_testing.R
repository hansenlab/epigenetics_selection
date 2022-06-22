library(Rfast)
getDM <- function(sorted_epiallele_vector){ #epiallele vector must be sorted based on the loeufs of the genes. 1st entry corresponds to the gene with the lowest loeuf etc
  pairs_matrix <- comb_n(sorted_epiallele_vector, 2)
  sum(pairs_matrix[1, ] - pairs_matrix[2, ])/dim(pairs_matrix)[2]
}

###proximal promoter
epiallele_vector <- meth_percentage
epiallele_vector[which(epiallele_vector >= 0.8)] <- 1
epiallele_vector[which(epiallele_vector <= 0.2)] <- 0
epiallele_vector <- epiallele_vector[which(epiallele_vector == 0 | epiallele_vector == 1)]
dm_proximal_prom <- getDM(epiallele_vector)

permutations_proximal_prom <- replicate(1000, {
  permuted_epiallele_vector <- epiallele_vector[sample(1:length(epiallele_vector), length(epiallele_vector), replace = FALSE)]
  getDM(permuted_epiallele_vector)
})

###promoter boundary
epiallele_vector <- hypometh_size
dm_hypometh_size <- getDM(epiallele_vector)

permutations_hypometh_size <- replicate(1000, {
  permuted_epiallele_vector <- epiallele_vector[sample(1:length(epiallele_vector), length(epiallele_vector), replace = FALSE)]
  getDM(permuted_epiallele_vector)
})


###tx end
epiallele_vector <- meth_percentage_tx_end
epiallele_vector[which(epiallele_vector >= 0.8)] <- 1
epiallele_vector[which(epiallele_vector <= 0.2)] <- 0
epiallele_vector <- epiallele_vector[which(epiallele_vector == 0 | epiallele_vector == 1)]
dm_tx_end <- getDM(epiallele_vector)

permutations_tx_end <- replicate(1000, {
  permuted_epiallele_vector <- epiallele_vector[sample(1:length(epiallele_vector), length(epiallele_vector), replace = FALSE)]
  getDM(permuted_epiallele_vector)
})

###h3k36me3
epiallele_vector <- df_k36me3_vs_expr$has_h3k36me3_peak
epiallele_vector[which(epiallele_vector == "yes")] <- 1
epiallele_vector[which(epiallele_vector == "no")] <- 0
epiallele_vector <- as.numeric(epiallele_vector) #not sure why I need to do that here
dm_h3k36me3 <- getDM(epiallele_vector)

permutations_h3k36me3 <- replicate(1000, {
  permuted_epiallele_vector <- epiallele_vector[sample(1:length(epiallele_vector), length(epiallele_vector), replace = FALSE)]
  getDM(permuted_epiallele_vector)
})

save(permutations_proximal_prom, permutations_hypometh_size, permutations_tx_end, permutations_tx_end_2, permutations_h3k36me3, 
     file = "permutations_selection_testing.rda") #tx_end_2 is using all coding transcripts, whereas tx_end is only for the transcripts for which we tested promoter meth also


###The following tests for h3k4me3 are done using the ucsd histone chip-seq data ("UCSD.H1.H3K4me3.LL312.narrowPeak.gz") from AnnotationHub
###h3k4me3 intensity (continuous)
epiallele_vector <- h3k4me3_intensity
dm_h3k4me3_intensity <- getDM(epiallele_vector)

permutations_h3k4me3_intensity <- replicate(1000, {
  permuted_epiallele_vector <- epiallele_vector[sample(1:length(epiallele_vector), length(epiallele_vector), replace = FALSE)]
  getDM(permuted_epiallele_vector)
})

###h3k4me3 presence (binary)
epiallele_vector <- proms_with_signal$signal_strength
epiallele_vector[which(epiallele_vector > 0)] <- 1
epiallele_vector[which(epiallele_vector == 0)] <- 0
epiallele_vector <- epiallele_vector[order(proms_with_signal$oe_lof_upper)]
dm_h3k4me3_presence <- getDM(epiallele_vector)

permutations_h3k4me3_presence <- replicate(1000, {
  permuted_epiallele_vector <- epiallele_vector[sample(1:length(epiallele_vector), length(epiallele_vector), replace = FALSE)]
  getDM(permuted_epiallele_vector)
})

save(permutations_h3k4me3_intensity, permutations_h3k4me3_presence, file = "permutations_h3k4me3.rda")



###figure
quartz(file = "permutations_proximal_prom.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
dm <- dm_proximal_prom
permutations <- permutations_proximal_prom
hist(permutations, breaks = 40, xlim = c(dm, max(permutations)), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(-0.07, -0.035, 0))
axis(2, at = c(0, 60, 120))
dev.off()

quartz(file = "permutations_hypometh_size.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
dm <- dm_hypometh_size
permutations <- permutations_hypometh_size
hist(permutations, breaks = 40, xlim = c(min(permutations), dm), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
#legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
#                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(0, 280, 560))
axis(2, at = c(0, 0.006, 0.012))
dev.off()

quartz(file = "permutations_tx_end.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
dm <- dm_tx_end
permutations <- permutations_tx_end
hist(permutations, breaks = 40, xlim = c(min(permutations), max(permutations)), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
#legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
#                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(- 0.012, 0, 0.012))
axis(2, at = c(0, 40, 80))
dev.off()

quartz(file = "permutations_h3k36me3.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
dm <- dm_h3k36me3
permutations <- permutations_h3k36me3
hist(permutations, breaks = 40, xlim = c(min(permutations), dm), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
#legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
#                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(0, 0.045, 0.09))
axis(2, at = c(0, 45, 90))
dev.off()

quartz(file = "permutations_h3k4me3_presence.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
dm <- dm_h3k4me3_presence
permutations <- permutations_h3k4me3_presence
hist(permutations, breaks = 40, xlim = c(min(permutations), dm), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
#legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
#                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(0, 0.06, 0.12))
axis(2, at = c(0, 50, 100))
dev.off()

quartz(file = "permutations_h3k4me3_intensity.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
dm <- dm_h3k4me3_intensity
permutations <- permutations_h3k4me3_intensity
hist(permutations, breaks = 40, xlim = c(min(permutations), dm), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
#legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
#                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(0, 2.25, 4.5))
axis(2, at = c(0, 0.8, 1.6))
dev.off()


####
epiallele_vector <- mean_meth
dm_mean_meth <- getDM(epiallele_vector)

#
epiallele_vector <- mean_meth_tx_end
dm_tx_end_mean <- getDM(epiallele_vector)
permutations_tx_end_mean <- replicate(1000, {
  permuted_epiallele_vector <- epiallele_vector[sample(1:length(epiallele_vector), length(epiallele_vector), replace = FALSE)]
  getDM(permuted_epiallele_vector)
})



