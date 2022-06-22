###
getAdjustedDM <- function(epiallele_vector, loeuf_vector, confounder_vector, epiallele_type = c("binary", "continuous")){
  model1 <- lm(loeuf_vector ~ confounder_vector)
  adjusted_loeuf <- residuals(model1)
  
  if (epiallele_type == "binary"){
    model2 <- glm(as.factor(epiallele_vector) ~ confounder_vector, family = "binomial")
    adjusted_epiallele <- residuals(model2, type = "deviance")
    adjusted_epiallele[which(adjusted_epiallele > 0)] <- 1
    adjusted_epiallele[which(adjusted_epiallele < 0)] <- 0
  } else if (epiallele_type == "continuous"){
    model2 <- lm(epiallele_vector ~ confounder_vector)
    adjusted_epiallele <- residuals(model2)
  }
  adjusted_epiallele <- adjusted_epiallele[order(adjusted_loeuf)]
  adjusted_dm <- getDM(adjusted_epiallele)
  adjusted_permutations <- replicate(1000, {
    permuted_epiallele_vector <- adjusted_epiallele[sample(1:length(adjusted_epiallele), length(adjusted_epiallele), replace = FALSE)]
    getDM(permuted_epiallele_vector)
  })
  list(adjusted_dm, adjusted_permutations)
}

getAdjustedDMWithoutPerm <- function(epiallele_vector, loeuf_vector, confounder_vector, epiallele_type = c("binary", "continuous")){
  model1 <- lm(loeuf_vector ~ confounder_vector)
  adjusted_loeuf <- residuals(model1)
  
  if (epiallele_type == "binary"){
    model2 <- glm(as.factor(epiallele_vector) ~ confounder_vector, family = "binomial")
    adjusted_epiallele <- residuals(model2, type = "deviance")
    adjusted_epiallele[which(adjusted_epiallele > 0)] <- 1
    adjusted_epiallele[which(adjusted_epiallele < 0)] <- 0
  } else if (epiallele_type == "continuous"){
    model2 <- lm(epiallele_vector ~ confounder_vector)
    adjusted_epiallele <- residuals(model2)
  }
  adjusted_epiallele <- adjusted_epiallele[order(adjusted_loeuf)]
  getDM(adjusted_epiallele)
}


###
###proximal promoter
#germline expr
epiallele_vector <- meth_percentage
epiallele_vector[which(epiallele_vector >= 0.8)] <- 1
epiallele_vector[which(epiallele_vector <= 0.2)] <- 0
epiallele_vector <- epiallele_vector[which(epiallele_vector == 0 | epiallele_vector == 1)]
loeufs <- human_loeufs[which(meth_percentage <= 0.2 | meth_percentage >= 0.8)]
gene_ids <- sapply(names(epiallele_vector), function(xx) proms_good_power$gene_id[which(proms_good_power$tx_id == xx)])
expr <- log2(median_expr_testis[gene_ids] + 1)
adjusted_dm_proximal_prom_germline <- getAdjustedDM(epiallele_vector, loeufs, expr, epiallele_type = "binary")

#esc expr
expr <- rowMeans(expr_h1_esc)[gene_ids]
adjusted_dm_proximal_prom_esc <- getAdjustedDM(epiallele_vector, loeufs, expr, epiallele_type = "binary")

#fetal expression
adjusted_dm_proximal_prom_fetal <- sapply(1:dim(fetal_expr)[2], function(xx) {
  expr <- fetal_expr[gene_ids, xx]
  getAdjustedDMWithoutPerm(epiallele_vector, loeufs, expr, epiallele_type = "binary")
})


###promoter boundary
#germline expr
epiallele_vector <- hypometh_size
loeufs <- tss_with_hypometh$oe_lof_upper[order(tss_with_hypometh$oe_lof_upper)]
gene_ids <- sapply(names(epiallele_vector), function(xx) proms_good_power$gene_id[which(proms_good_power$tx_id == xx)])
expr <- log2(median_expr_testis[gene_ids] + 1)
adjusted_dm_hypometh_size_germline <- getAdjustedDM(epiallele_vector, loeufs, expr, epiallele_type = "continuous")

#esc expr
expr <- rowMeans(expr_h1_esc)[gene_ids]
adjusted_dm_hypometh_size_esc <- getAdjustedDM(epiallele_vector, loeufs, expr, epiallele_type = "continuous")

#fetal expression
adjusted_dm_hypometh_size_fetal <- sapply(1:dim(fetal_expr)[2], function(xx) {
  expr <- fetal_expr[gene_ids, xx]
  getAdjustedDMWithoutPerm(epiallele_vector, loeufs, expr, epiallele_type = "continuous")
})

###H3K36me3
#germline expr
epiallele_vector <- df_k36me3_vs_expr$has_h3k36me3_peak
epiallele_vector[which(epiallele_vector == "yes")] <- 1
epiallele_vector[which(epiallele_vector == "no")] <- 0
epiallele_vector <- as.numeric(epiallele_vector) #not sure why I need to do that here
loeufs <- df_k36me3_vs_expr$loeuf
gene_ids <- sapply(rownames(df_k36me3_vs_expr), function(xx) tx_all$gene_id[which(tx_all$tx_id == xx)])
expr <- log2(median_expr_testis[gene_ids] + 1)
adjusted_dm_h3k36me3_germline <- getAdjustedDM(epiallele_vector, loeufs, expr, epiallele_type = "binary")

#esc expr
expr <- rowMeans(expr_h1_esc)[gene_ids]
adjusted_dm_h3k36me3_esc <- getAdjustedDM(epiallele_vector, loeufs, expr, epiallele_type = "binary")

#fetal expression
adjusted_dm_h3k36me3_fetal <- sapply(1:dim(fetal_expr)[2], function(xx) {
  expr <- fetal_expr[gene_ids, xx]
  getAdjustedDMWithoutPerm(epiallele_vector, loeufs, expr, epiallele_type = "binary")
})

###H3K4me3 presence
#esc expr (no need to check germline expr here since the chip seq data are from H1 embryonic stem cells)
epiallele_vector <- proms_with_signal$signal_strength
epiallele_vector[which(epiallele_vector > 0)] <- 1
epiallele_vector[which(epiallele_vector == 0)] <- 0
epiallele_vector <- epiallele_vector[order(proms_with_signal$oe_lof_upper)]
loeufs <- proms_with_signal$oe_lof_upper[order(proms_with_signal$oe_lof_upper)]
gene_ids <- proms_with_signal$gene_id[order(proms_with_signal$oe_lof_upper)]
expr <- rowMeans(expr_h1_esc)[gene_ids]

adjusted_dm_h3k4me3_presence <- getAdjustedDM(epiallele_vector, loeufs, expr, epiallele_type = "binary")

#fetal expression
adjusted_dm_h3k4me3_presence_fetal <- sapply(1:dim(fetal_expr)[2], function(xx) {
  expr <- fetal_expr[gene_ids, xx]
  getAdjustedDMWithoutPerm(epiallele_vector, loeufs, expr, epiallele_type = "binary")
})

###H3K4me3 signal
#esc expr
epiallele_vector <- h3k4me3_intensity
loeufs <- proms_with_signal_2$oe_lof_upper[order(proms_with_signal_2$oe_lof_upper)]
gene_ids <- proms_with_signal_2$gene_id[order(proms_with_signal_2$oe_lof_upper)]
expr <- rowMeans(expr_h1_esc)[gene_ids]
adjusted_dm_h3k4me3_intensity <- getAdjustedDM(epiallele_vector, loeufs, expr, epiallele_type = "continuous")

#fetal expression
adjusted_dm_h3k4me3_intensity_fetal <- sapply(1:dim(fetal_expr)[2], function(xx) {
  expr <- fetal_expr[gene_ids, xx]
  getAdjustedDMWithoutPerm(epiallele_vector, loeufs, expr, epiallele_type = "continuous")
})


save(adjusted_dm_proximal_prom_germline, adjusted_dm_hypometh_size_germline, adjusted_dm_h3k36me3_germline, 
     adjusted_dm_proximal_prom_esc, adjusted_dm_hypometh_size_esc, adjusted_dm_h3k36me3_esc, adjusted_dm_h3k4me3_presence, 
     adjusted_dm_h3k4me3_intensity, file = "adjusted_epiallele_testing.rda")

save(adjusted_dm_proximal_prom_fetal, adjusted_dm_hypometh_size_fetal, adjusted_dm_h3k36me3_fetal, adjusted_dm_h3k4me3_presence_fetal, 
     adjusted_dm_h3k4me3_intensity_fetal, file = "adjusted_epiallele_testing_fetal_expr.rda")

quartz(file = "permutations_adjusted_epialleles_germline.pdf", height = 2.2, width = 6.5, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
dm <- adjusted_dm_proximal_prom_germline[[1]]
permutations <- adjusted_dm_proximal_prom_germline[[2]]
hist(permutations, breaks = 40, xlim = c(dm, max(permutations)), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "adjusted dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(-0.05, -0.025, 0))
axis(2, at = c(0, 70, 140))

dm <- adjusted_dm_hypometh_size_germline[[1]]
permutations <- adjusted_dm_hypometh_size_germline[[2]]
hist(permutations, breaks = 40, xlim = c(min(permutations), dm), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "adjusted dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
#legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
#                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(0, 300, 600))
axis(2, at = c(0, 0.007, 0.014))

dm <- adjusted_dm_h3k36me3_germline[[1]]
permutations <- adjusted_dm_h3k36me3_germline[[2]]
hist(permutations, breaks = 40, xlim = c(min(permutations), dm), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "adjusted dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
#legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
#                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(0, 0.025, 0.05))
axis(2, at = c(0, 75, 150))
dev.off()

quartz(file = "permutations_adjusted_epialleles_esc.pdf", height = 2.2, width = 11, pointsize = 8, type = "pdf")
par(mfrow = c(1, 5))
dm <- adjusted_dm_proximal_prom_esc[[1]]
permutations <- adjusted_dm_proximal_prom_esc[[2]]
hist(permutations, breaks = 40, xlim = c(dm, max(permutations)), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "adjusted dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(-0.03, -0.015, 0))
axis(2, at = c(0, 75, 150))

dm <- adjusted_dm_hypometh_size_esc[[1]]
permutations <- adjusted_dm_hypometh_size_esc[[2]]
hist(permutations, breaks = 40, xlim = c(min(permutations), dm), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "adjusted dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
#legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
#                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(0, 300, 600))
axis(2, at = c(0, 0.007, 0.014))

dm <- adjusted_dm_h3k36me3_esc[[1]]
permutations <- adjusted_dm_h3k36me3_esc[[2]]
hist(permutations, breaks = 40, xlim = c(min(permutations), dm), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "adjusted dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
#legend <- legend("top", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
#                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(1, at = c(0, 0.025, 0.05))
axis(2, at = c(0, 60, 120))

dm <- adjusted_dm_h3k4me3_presence[[1]]
permutations <- adjusted_dm_h3k4me3_presence[[2]]
hist(permutations, breaks = 40, xlim = c(min(permutations), dm), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "adjusted dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
axis(1, at = c(0, 0.015, 0.03))
axis(2, at = c(0, 40, 80))

dm <- adjusted_dm_h3k4me3_intensity[[1]]
permutations <- adjusted_dm_h3k4me3_intensity[[2]]
hist(permutations, breaks = 40, xlim = c(min(permutations), dm), lty = 0,
     main = "", font.main = 1, freq = FALSE, 
     col = "cornflowerblue", xlab = "adjusted dM", xaxt = 'n', yaxt = 'n')
abline(v = dm, lwd = 3, col = alpha("red", 0.75))
axis(1, at = c(0, 1, 2))
axis(2, at = c(0, 0.75, 1.5))
dev.off()

quartz(file = "dM_adjusted_epialleles_fetal_expression.pdf", height = 2.2, width = 11, pointsize = 8, type = "pdf")
par(mfrow = c(1, 5))
plot(adjusted_dm_proximal_prom_fetal, pch = 19, cex = 0.8, col = alpha("red", 0.75), bty = 'l', 
     ylim = c(min(adjusted_dm_proximal_prom_fetal), 0.01),
     main = "DNA methylation (proximal prom)", font.main = 1, ylab = "dM adjusted for expression", xlab = "172 fetal tissues", 
     xaxt = 'n', yaxt = 'n')
axis(2, at = c(-0.07, 0))

plot(adjusted_dm_hypometh_size_fetal, pch = 19, cex = 0.8, col = alpha("red", 0.75), bty = 'l', 
     ylim = c(-0.01, max(adjusted_dm_hypometh_size_fetal)),
     main = "hypomethylated region size", font.main = 1, ylab = "dM adjusted for expression", xlab = "172 fetal tissues", 
     xaxt = 'n', yaxt = 'n')
axis(2, at = c(0, 550))

plot(adjusted_dm_h3k36me3_fetal, pch = 19, cex = 0.8, col = alpha("red", 0.75), bty = 'l', 
     ylim = c(-0.01, max(adjusted_dm_h3k36me3_fetal)),
     main = "H3K36me3 presence (coding seq)", font.main = 1, ylab = "dM adjusted for expression", xlab = "172 fetal tissues", 
     xaxt = 'n', yaxt = 'n')
axis(2, at = c(0, 0.08))

plot(adjusted_dm_h3k4me3_presence_fetal, pch = 19, cex = 0.8, col = alpha("red", 0.75), bty = 'l', 
     ylim = c(-0.01, max(adjusted_dm_h3k4me3_presence_fetal)),
     main = "H3K4me3 presence (proximal prom)", font.main = 1, ylab = "dM adjusted for expression", xlab = "172 fetal tissues", 
     xaxt = 'n', yaxt = 'n')
axis(2, at = c(0, 0.1))

plot(adjusted_dm_h3k4me3_intensity_fetal, pch = 19, cex = 0.8, col = alpha("red", 0.75), bty = 'l', 
     ylim = c(-0.01, max(adjusted_dm_h3k4me3_intensity_fetal)),
     main = "H3K4me3 signal intensity\n(proximal prom)", font.main = 1, ylab = "dM adjusted for expression", xlab = "172 fetal tissues", 
     xaxt = 'n', yaxt = 'n')
axis(2, at = c(0, 4))
dev.off()



###boxplots instead
quartz(file = "permutations_adjusted_epialleles_proximal_prom_1.pdf", height = 2.2, width = 1.1, pointsize = 8, type = "pdf")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
dm1 <- adjusted_dm_proximal_prom_germline[[1]]
permutations1 <- adjusted_dm_proximal_prom_germline[[2]]

boxplot(permutations1, frame = FALSE, main = "proximal prom methylation", font.main = 1, cex.main = 0.57,
        lty = "solid", col = "cornflowerblue", yaxt = 'n', ylim = c(dm1, max(permutations1)),
        xlim = c(0.8, 2.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "adjusted dM")
points(1, dm1, cex = 1.25, col = alpha("red", 0.75), pch = 19)
axis(1, at = 1, labels = c("adjusted for\ngermline expr"))
axis(2, at = c(-0.05, -0.025, 0), las = 2)
dev.off()

quartz(file = "permutations_adjusted_epialleles_proximal_prom_2.pdf", height = 2.2, width = 1.1, pointsize = 8, type = "pdf")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
dm2 <- adjusted_dm_proximal_prom_esc[[1]]
permutations2 <- adjusted_dm_proximal_prom_esc[[2]]

boxplot(permutations2, frame = FALSE, main = "proximal prom methylation", font.main = 1, cex.main = 0.57,
        lty = "solid", col = "cornflowerblue", yaxt = 'n', ylim = c(dm1, max(permutations2)),
        xlim = c(0.8, 2.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "adjusted dM")
points(1, dm2, cex = 1.25, col = alpha("red", 0.75), pch = 19)
#legend <- legend("bottom", legend = c("null (permutations)", "observed"), bty = 'n', cex = 0.75, 
#                 lwd = 2.5, col = c("cornflowerblue", alpha("red", 0.75)))
axis(2, at = c(-0.05, -0.025, 0), las = 2)
axis(1, at = 1, labels = c("adjusted for\nESC expr"))
dev.off()

quartz(file = "permutations_adjusted_epialleles_hypometh_size_1.pdf", height = 2.2, width = 1.1, pointsize = 8, type = "pdf")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
dm1 <- adjusted_dm_hypometh_size_germline[[1]]
permutations1 <- adjusted_dm_hypometh_size_germline[[2]]

boxplot(permutations1, frame = FALSE, main = "size of hypometh region", font.main = 1, cex.main = 0.57,
        lty = "solid", col = "cornflowerblue", yaxt = 'n', ylim = c(min(permutations1), dm1),
        xlim = c(0.8, 2.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "adjusted dM")
points(1, dm1, cex = 1.25, col = alpha("red", 0.75), pch = 19)
axis(2, at = c(0, 300, 600), las = 2)
axis(1, at = 1, labels = c("adjusted for\ngermline expr"))
dev.off()

quartz(file = "permutations_adjusted_epialleles_hypometh_size_2.pdf", height = 2.2, width = 1.1, pointsize = 8, type = "pdf")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
dm2 <- adjusted_dm_hypometh_size_esc[[1]]
permutations2 <- adjusted_dm_hypometh_size_esc[[2]]
boxplot(permutations2, frame = FALSE, main = "size of hypometh region", font.main = 1, cex.main = 0.57,
        lty = "solid", col = "cornflowerblue", yaxt = 'n', ylim = c(min(permutations2), dm2),
        xlim = c(0.8, 2.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "adjusted dM")
points(1, dm2, cex = 1.25, col = alpha("red", 0.75), pch = 19)
axis(2, at = c(0, 300, 600), las = 2)
axis(1, at = 1, labels = c("adjusted for\nESC expr"))
dev.off()


quartz(file = "permutations_adjusted_epialleles_h3k36me3_1.pdf", height = 2.2, width = 1.1, pointsize = 8, type = "pdf")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
dm1 <- adjusted_dm_h3k36me3_germline[[1]]
permutations1 <- adjusted_dm_h3k36me3_germline[[2]]

boxplot(permutations1, frame = FALSE, main = "coding H3K36me3 peak", font.main = 1, cex.main = 0.57,
        lty = "solid", col = "cornflowerblue", yaxt = 'n', ylim = c(min(permutations1), dm1),
        xlim = c(0.8, 2.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "adjusted dM")
points(1, dm1, cex = 1.25, col = alpha("red", 0.75), pch = 19)
axis(2, at = c(0, 0.025, 0.05), las = 2)
axis(1, at = 1, labels = c("adjusted for\ngermline expr"))
dev.off()

quartz(file = "permutations_adjusted_epialleles_h3k36me3_2.pdf", height = 2.2, width = 1.1, pointsize = 8, type = "pdf")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
dm2 <- adjusted_dm_h3k36me3_esc[[1]]
permutations2 <- adjusted_dm_h3k36me3_esc[[2]]
boxplot(permutations2, frame = FALSE, main = "coding H3K36me3 peak", font.main = 1, cex.main = 0.57,
        lty = "solid", col = "cornflowerblue", yaxt = 'n', ylim = c(min(permutations2), dm2),
        xlim = c(0.8, 2.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "adjusted dM")
points(1, dm2, cex = 1.25, col = alpha("red", 0.75), pch = 19)
axis(2, at = c(0, 0.025, 0.05), las = 2)
axis(1, at = 1, labels = c("adjusted for\nESC expr"))
dev.off()

quartz(file = "permutations_adjusted_epialleles_h3k4me3_presence.pdf", height = 2.2, width = 1.1, pointsize = 8, type = "pdf")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
dm1 <- adjusted_dm_h3k4me3_presence[[1]]
permutations1 <- adjusted_dm_h3k4me3_presence[[2]]

boxplot(permutations1, frame = FALSE, main = "proximal prom\nH3K4me3 peak", font.main = 1, cex.main = 0.57,
        lty = "solid", col = "cornflowerblue", yaxt = 'n', ylim = c(min(permutations1), dm1),
        xlim = c(0.8, 2.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "adjusted dM")
points(1, dm1, cex = 1.25, col = alpha("red", 0.75), pch = 19)
axis(2, at = c(0, 0.015, 0.03), las = 2)
axis(1, at = 1, labels = c("adjusted for\nESC expr"))
dev.off()

quartz(file = "permutations_adjusted_epialleles_h3k4me3_intensity.pdf", height = 2.2, width = 1.1, pointsize = 8, type = "pdf")
par(font.main = 1, mar = c(3,4,1,1) + 0.1)
dm2 <- adjusted_dm_h3k4me3_intensity[[1]]
permutations2 <- adjusted_dm_h3k4me3_intensity[[2]]
boxplot(permutations2, frame = FALSE, main = "proximal prom\nH3K4me3 intensity", font.main = 1, cex.main = 0.57,
        lty = "solid", col = "cornflowerblue", yaxt = 'n', ylim = c(min(permutations2), dm2),
        xlim = c(0.8, 2.2), medlty = 1, medlwd = 0.8, at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "adjusted dM")
points(1, dm2, cex = 1.25, col = alpha("red", 0.75), pch = 19)
axis(2, at = c(0, 1, 2), las = 2)
axis(1, at = 1, labels = c("adjusted for\nESC expr"))
dev.off()



