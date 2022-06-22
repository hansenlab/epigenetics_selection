addSignalToPromGrange <- function(proms_grange, signal_grange){
  overlaps <- findOverlaps(proms_grange, signal_grange)
  proms_grange$signal_strength <- 0
  peak_by_prom <- split(signal_grange$signalValue[subjectHits(overlaps)], queryHits(overlaps)) #ensure there is a column named signalValue
  proms_grange$signal_strength[unique(queryHits(overlaps))] <- sapply(peak_by_prom, mean) 
  proms_grange
}


ah <- AnnotationHub()
ah <- subset(ah, species == "Homo sapiens")
qhs <- query(ah, "H3K4me3")

qhs$title[grep("narrow", qhs$title)]
#gr1 <- subset(qhs, title == "BI.H1.H3K4me3.Solexa-8038.narrowPeak.gz")[[1]]
gr1 <- subset(qhs, title == "UCSD.H1.H3K4me3.LL312.narrowPeak.gz")[[1]]
#gr2 <- subset(qhs, title == "UCSD.H1_Derived_Neuronal_Progenitor_Cultured_Cells.H3K4me3.AK199.narrowPeak.gz")[[1]]

proms_with_signal <- addSignalToPromGrange(proms_good_power-1500, gr1)
proms_with_signal_2 <- proms_with_signal[which(proms_with_signal$signal_strength > 0)]

#######
proms_granges_for_plot <- proms_with_signal_2
proms_granges_for_plot$group <- 10
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.9) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.8))] <- 9
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.8) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.7))] <- 8
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.7) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.6))] <- 7
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.6) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.5))] <- 6
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.5) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.4))] <- 5
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.4) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.3))] <- 4
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.3) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.2))] <- 3
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.2) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.1))] <- 2
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.1) 
)] <- 1

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))


proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                signal = proms_granges_for_plot$signal_strength)


quartz(file = "h3k4me3_signal_value_loeuf_ridgeplot.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
values_high_decile <- proms_df_for_plot$signal[which(proms_df_for_plot$group1 == 10)]
mode_to_plot <- density(values_high_decile)$x[which.max(density(values_high_decile)$y)]
ggplot(proms_df_for_plot, aes(x = signal, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), scale = 1.1, from = 0, to = 80) + 
  labs(x = "h3k4me3 signal value", y = "LOEUF decile") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size = 8)) + theme(axis.title=element_text(size = 10)) + geom_vline(xintercept = mode_to_plot)
  #geom_vline(xintercept = median(proms_df_for_plot$signal[which(proms_df_for_plot$group1 == 10)]))
dev.off()

quartz(file = "h3k4me3_signal_value_loeuf.pdf", height = 5, width = 3, type = "pdf")
ggplot(proms_df_for_plot, aes(x = group1, y = signal)) + 
  geom_boxplot(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(1, 20)) + labs(y = "h3k4me3 signal value", x = "LOEUF decile") + 
  theme(axis.text=element_text(size=10)) + theme_classic()
dev.off()

###
quartz(file = "h3k4me3_vs_loeuf_promoter.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(1, getPercentOfPromotersWithSignal(proms_with_signal, proms_with_signal$oe_lof_upper, proms_good_power$oe_lof_upper, 0.1, "less"), pch = 19, col = alpha("red", 0.59), 
     bty = 'l', xlab = "LOEUF decile", cex = 1.2, ylab = "% genes w/ H3K4me3 peak", main = "proximal promoter\n(+/- 500bp from TSS)", 
     xlim = c(0.8, 10.2), ylim = c(0.6, 1), yaxt = 'n', xaxt= 'n', font.main = 1)
points(2, getPercentOfPromotersWithSignal(proms_with_signal, proms_with_signal$oe_lof_upper, proms_with_signal$oe_lof_upper, c(0.1, 0.2), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(3, getPercentOfPromotersWithSignal(proms_with_signal, proms_with_signal$oe_lof_upper, proms_with_signal$oe_lof_upper, c(0.2, 0.3), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(4, getPercentOfPromotersWithSignal(proms_with_signal, proms_with_signal$oe_lof_upper, proms_with_signal$oe_lof_upper, c(0.3, 0.4), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(5, getPercentOfPromotersWithSignal(proms_with_signal, proms_with_signal$oe_lof_upper, proms_with_signal$oe_lof_upper, c(0.4, 0.5), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(6, getPercentOfPromotersWithSignal(proms_with_signal, proms_with_signal$oe_lof_upper, proms_with_signal$oe_lof_upper, c(0.5, 0.6), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(7, getPercentOfPromotersWithSignal(proms_with_signal, proms_with_signal$oe_lof_upper, proms_with_signal$oe_lof_upper, c(0.6, 0.7), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(8, getPercentOfPromotersWithSignal(proms_with_signal, proms_with_signal$oe_lof_upper, proms_with_signal$oe_lof_upper, c(0.7, 0.8), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(9, getPercentOfPromotersWithSignal(proms_with_signal, proms_with_signal$oe_lof_upper, proms_with_signal$oe_lof_upper, c(0.8, 0.9), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
points(10, getPercentOfPromotersWithSignal(proms_with_signal, proms_with_signal$oe_lof_upper, proms_with_signal$oe_lof_upper, 0.9, "greater"), pch = 19, col = alpha("red", 0.59), cex = 1.2)
axis(1, at = 1:10, cex.axis = 0.7)
axis(2, at = c(0.6, 1))
dev.off()




###relationship with expression
proms_with_signal_2$expr_esc <- rowMeans(expr_h1_esc[proms_with_signal_2$gene_id, ])

proms_granges_for_plot <- proms_with_signal_2
proms_granges_for_plot$group <- 10
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_esc < quantile(proms_granges_for_plot$expr_esc, 0.9) 
                                   & proms_granges_for_plot$expr_esc > quantile(proms_granges_for_plot$expr_esc, 0.8))] <- 9
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_esc < quantile(proms_granges_for_plot$expr_esc, 0.8) 
                                   & proms_granges_for_plot$expr_esc > quantile(proms_granges_for_plot$expr_esc, 0.7))] <- 8
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_esc < quantile(proms_granges_for_plot$expr_esc, 0.7) 
                                   & proms_granges_for_plot$expr_esc > quantile(proms_granges_for_plot$expr_esc, 0.6))] <- 7
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_esc < quantile(proms_granges_for_plot$expr_esc, 0.6) 
                                   & proms_granges_for_plot$expr_esc > quantile(proms_granges_for_plot$expr_esc, 0.5))] <- 6
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_esc < quantile(proms_granges_for_plot$expr_esc, 0.5) 
                                   & proms_granges_for_plot$expr_esc > quantile(proms_granges_for_plot$expr_esc, 0.4))] <- 5
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_esc < quantile(proms_granges_for_plot$expr_esc, 0.4) 
                                   & proms_granges_for_plot$expr_esc > quantile(proms_granges_for_plot$expr_esc, 0.3))] <- 4
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_esc < quantile(proms_granges_for_plot$expr_esc, 0.3) 
                                   & proms_granges_for_plot$expr_esc > quantile(proms_granges_for_plot$expr_esc, 0.2))] <- 3
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_esc < quantile(proms_granges_for_plot$expr_esc, 0.2) 
                                   & proms_granges_for_plot$expr_esc > quantile(proms_granges_for_plot$expr_esc, 0.1))] <- 2
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_esc < quantile(proms_granges_for_plot$expr_esc, 0.1) 
)] <- 1

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))


proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                signal = proms_granges_for_plot$signal_strength)


quartz(file = "h3k4me3_signal_value_expr_esc.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
values_high_decile <- proms_df_for_plot$signal[which(proms_df_for_plot$group1 == 10)]
mode_to_plot <- density(values_high_decile)$x[which.max(density(values_high_decile)$y)]
ggplot(proms_df_for_plot, aes(x = signal, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), scale = 1.1, from = 0, to = 80) + 
  labs(x = "h3k4me3 signal value", y = "ESC expression decile") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size = 8)) + theme(axis.title=element_text(size = 10)) + geom_vline(xintercept = mode_to_plot)
#geom_vline(xintercept = median(proms_df_for_plot$signal[which(proms_df_for_plot$group1 == 10)]))
dev.off()



##
h3k4me3_intensity <- proms_with_signal_2$signal_strength
h3k4me3_intensity <- h3k4me3_intensity[order(proms_with_signal_2$oe_lof_upper)]

####
quartz(file = "h3k4me3_intensity_histogram.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
hist(proms_with_signal_2$signal_strength, col = alpha("red", 0.6), lty = 0, 
     breaks = 65, freq = FALSE, xlab = "H3K4me3 signal value", cex.lab = 1.2, 
     yaxt = 'n', xaxt = 'n', main = "", font.main = 1, cex.main = 1.25)
axis(1, at = c(0, 50, 100), cex.axis = 1.2)
axis(2, at = c(0, 0.025), cex.axis = 1.2)
#abline(v = c(40, 80), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()





