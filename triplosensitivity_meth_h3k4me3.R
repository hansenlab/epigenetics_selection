triplosensitivity <- read_csv('haplo_triplo.csv')
triplosensitivity <- triplosensitivity[which(triplosensitivity$gene %in% constraint$gene[which(constraint$canonical == TRUE)]), ]
triplosensitivity$tx_id <- sapply(triplosensitivity$gene, 
                                  function(xx) constraint$transcript[which(constraint$gene == xx & constraint$canonical == TRUE)])
triplosensitivity <- triplosensitivity[which(lengths(triplosensitivity$tx_id) == 1), ]
triplosensitivity$tx_id <- unlist(triplosensitivity$tx_id)
triplosensitivity <- as.data.frame(triplosensitivity)
rownames(triplosensitivity) <- triplosensitivity$tx_id
triplosensitivity1 <- triplosensitivity[names(meth_percentage), ]
#triplosensitivity1 <- triplosensitivity1[which(triplosensitivity1$pHI <= 0.4), ]

meth_vs_loeuf_vec2 <- c(getPercentOfMethPromoters(meth_percentage, triplosensitivity1$pTS, 0.1, "less", 0.8), 
                        getPercentOfMethPromoters(meth_percentage, triplosensitivity1$pTS, c(0.1, 0.2), "in between", 0.8),
                        getPercentOfMethPromoters(meth_percentage, triplosensitivity1$pTS, c(0.2, 0.3), "in between", 0.8),
                        getPercentOfMethPromoters(meth_percentage, triplosensitivity1$pTS, c(0.3, 0.4), "in between", 0.8),
                        getPercentOfMethPromoters(meth_percentage, triplosensitivity1$pTS, c(0.4, 0.5), "in between", 0.8),
                        getPercentOfMethPromoters(meth_percentage, triplosensitivity1$pTS, c(0.5, 0.6), "in between", 0.8),
                        getPercentOfMethPromoters(meth_percentage, triplosensitivity1$pTS, c(0.6, 0.7), "in between", 0.8),
                        getPercentOfMethPromoters(meth_percentage, triplosensitivity1$pTS, c(0.7, 0.8), "in between", 0.8),
                        getPercentOfMethPromoters(meth_percentage, triplosensitivity1$pTS, c(0.8, 0.9), "in between", 0.8),
                        getPercentOfMethPromoters(meth_percentage, triplosensitivity1$pTS, 0.9, "greater", 0.8))

quartz(file = "methylation_vs_triplosensitivity_promoter.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(meth_vs_loeuf_vec2, pch = 19, cex = 1.2, col = alpha("red", 0.75), bty = 'l', ylim = c(0, 0.04),
     main = "proximal promoter\n(+/- 500bp from TSS)", font.main = 1, ylab = "% genes w/ methylation", xlab = "triplosensitivity decile", 
     xaxt = 'n', yaxt = 'n')
axis(1, at = 1:10, cex.axis = 0.7)
axis(2, at = c(0, 0.04))
dev.off()


triplosensitivity2 <- triplosensitivity[proms_with_signal$tx_id, ]
#triplosensitivity2 <- triplosensitivity2[which(triplosensitivity2$pHI <= 0.4), ]
quartz(file = "h3k4me3_vs_triplosensitivity_promoter.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(1, getPercentOfPromotersWithSignal(proms_with_signal, triplosensitivity2$pTS, triplosensitivity2$pTS, 0.1, "less"), pch = 19, col = alpha("red", 0.59), 
     bty = 'l', xlab = "triplosensitivity decile", cex = 1.2, ylab = "% genes w/ H3K4me3 peak", main = "proximal promoter\n(+/- 500bp from TSS)", 
     xlim = c(0.8, 10.2), ylim = c(0.5, 1), yaxt = 'n', xaxt= 'n', font.main = 1)
points(2, getPercentOfPromotersWithSignal(proms_with_signal, triplosensitivity2$pTS, triplosensitivity2$pTS, c(0.1, 0.2), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(3, getPercentOfPromotersWithSignal(proms_with_signal, triplosensitivity2$pTS, triplosensitivity2$pTS, c(0.2, 0.3), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(4, getPercentOfPromotersWithSignal(proms_with_signal, triplosensitivity2$pTS, triplosensitivity2$pTS, c(0.3, 0.4), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(5, getPercentOfPromotersWithSignal(proms_with_signal, triplosensitivity2$pTS, triplosensitivity2$pTS, c(0.4, 0.5), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(6, getPercentOfPromotersWithSignal(proms_with_signal, triplosensitivity2$pTS, triplosensitivity2$pTS, c(0.5, 0.6), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(7, getPercentOfPromotersWithSignal(proms_with_signal, triplosensitivity2$pTS, triplosensitivity2$pTS, c(0.6, 0.7), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(8, getPercentOfPromotersWithSignal(proms_with_signal, triplosensitivity2$pTS, triplosensitivity2$pTS, c(0.7, 0.8), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(9, getPercentOfPromotersWithSignal(proms_with_signal, triplosensitivity2$pTS, triplosensitivity2$pTS, c(0.8, 0.9), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(10, getPercentOfPromotersWithSignal(proms_with_signal, triplosensitivity2$pTS, triplosensitivity2$pTS, 0.9, "greater"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
axis(1, at = 1:10, cex.axis = 0.7)
axis(2, at = c(0.6, 1))
dev.off()

