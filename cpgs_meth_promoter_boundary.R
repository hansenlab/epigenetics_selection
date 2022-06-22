####same analysis of methylation vs LOEUF but using size of hypometh regions
edb <- EnsDb.Hsapiens.v75
tx_all <- transcripts(edb, columns = c("tx_id", "gene_id"))
tx_all <- tx_all[proms_good_power$tx_id]
tss_all <- resize(tx_all, 1, fix = "start")
tss_all$oe_lof_upper <- proms_good_power$oe_lof_upper
seqlevelsStyle(tss_all) <- "ucsc"
genome(tss_all) <- "hg19"


hypometh_regions_human <- import("~/Downloads/germline_methylation/human_germline_hypometh_regions.bed", 
                               format = "BED")
genome(hypometh_regions_human) <- "hg19"
hypometh_regions_human <- hypometh_regions_human[which(seqnames(hypometh_regions_human) %in% names(Hsapiens)[1:22])]
seqlevels(hypometh_regions_human) <- seqlevels(hypometh_regions_human)[1:22]

tss_all$hypometh_size_human <- 0
overlaps <- findOverlaps(tss_all, hypometh_regions_human)
tss_all$hypometh_size_human[unique(queryHits(overlaps))] <- width(hypometh_regions_human[subjectHits(overlaps)])

tss_all <- tss_all[proms_good_power_no_bidir$tx_id] #we exclude bidirectional promoters because these will artificially appear as large hypomethylated regions
tss_with_hypometh <- tss_all[which(tss_all$hypometh_size_human > 0)]


hypometh_size <- tss_with_hypometh$hypometh_size_human
names(hypometh_size) <- tss_with_hypometh$tx_id
hypometh_size <- hypometh_size[order(tss_with_hypometh$oe_lof_upper)]

######figure showing association with hypometh region size and loeuf
proms_granges_for_plot <- tss_with_hypometh
proms_granges_for_plot$group <- 10
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.9, na.rm = TRUE) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.8, na.rm = TRUE))] <- 9
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.8, na.rm = TRUE) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.7, na.rm = TRUE))] <- 8
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.7, na.rm = TRUE) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.6, na.rm = TRUE))] <- 7
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.6, na.rm = TRUE) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.5, na.rm = TRUE))] <- 6
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.5, na.rm = TRUE) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.4, na.rm = TRUE))] <- 5
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.4, na.rm = TRUE)
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.3, na.rm = TRUE))] <- 4
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.3, na.rm = TRUE) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.2, na.rm = TRUE))] <- 3
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.2, na.rm = TRUE) 
                                   & proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.1, na.rm = TRUE))] <- 2
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.1, na.rm = TRUE)) 
] <- 1

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))



proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                hypometh_size = proms_granges_for_plot$hypometh_size_human)

quartz(file = "hypometh_size_loeuf_ridgeplot.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
#values_high_decile <- proms_df_for_plot$hypometh_size[which(proms_df_for_plot$group1 == 10)]
#mode_to_plot <- density(values_high_decile)$x[which.max(density(values_high_decile)$y)]
ggplot(proms_df_for_plot, aes(x = hypometh_size, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), scale = 1.1, from = 0, to = 7500) + 
  labs(x = "hypometh region size", y = "LOEUF decile") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size = 8)) + theme(axis.title=element_text(size = 10)) +
  #geom_vline(xintercept = mode_to_plot)
  geom_vline(xintercept = median(proms_df_for_plot$hypometh_size[which(proms_df_for_plot$group1 == 10)]))
dev.off()


#####
quartz(file = "hypometh_size_loeuf.pdf", height = 2.4, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = group1, y = hypometh_size)) + 
  geom_boxplot(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(500, 8000)) + labs(y = "TSS hypometh region size", x = "LOEUF decile") + 
  theme(axis.text=element_text(size=10)) + theme_classic()
dev.off()




###explicitly visualize relationship with expression
#esc expr
tss_with_hypometh$expr_esc <- rowMeans(expr_h1_esc[tss_with_hypometh$gene_id, ])
proms_granges_for_plot <- tss_with_hypometh
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
                                hypometh_size = proms_granges_for_plot$hypometh_size_human)

quartz(file = "hypometh_size_expr_esc.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = hypometh_size, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), scale = 1.1, from = 0, to = 7500) + 
  labs(x = "hypometh region size", y = "ESC expression decile") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size = 8)) + theme(axis.title=element_text(size = 10)) +
  #geom_vline(xintercept = mode_to_plot)
  geom_vline(xintercept = median(proms_df_for_plot$hypometh_size[which(proms_df_for_plot$group1 == 10)]))

dev.off()



#germline expr
tss_with_hypometh$expr_germline <- log2(median_expr_testis[tss_with_hypometh$gene_id] + 1)
proms_granges_for_plot <- tss_with_hypometh
proms_granges_for_plot$group <- 10
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_germline < quantile(proms_granges_for_plot$expr_germline, 0.9) 
                                   & proms_granges_for_plot$expr_germline > quantile(proms_granges_for_plot$expr_germline, 0.8))] <- 9
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_germline < quantile(proms_granges_for_plot$expr_germline, 0.8) 
                                   & proms_granges_for_plot$expr_germline > quantile(proms_granges_for_plot$expr_germline, 0.7))] <- 8
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_germline < quantile(proms_granges_for_plot$expr_germline, 0.7) 
                                   & proms_granges_for_plot$expr_germline > quantile(proms_granges_for_plot$expr_germline, 0.6))] <- 7
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_germline < quantile(proms_granges_for_plot$expr_germline, 0.6) 
                                   & proms_granges_for_plot$expr_germline > quantile(proms_granges_for_plot$expr_germline, 0.5))] <- 6
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_germline < quantile(proms_granges_for_plot$expr_germline, 0.5) 
                                   & proms_granges_for_plot$expr_germline > quantile(proms_granges_for_plot$expr_germline, 0.4))] <- 5
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_germline < quantile(proms_granges_for_plot$expr_germline, 0.4) 
                                   & proms_granges_for_plot$expr_germline > quantile(proms_granges_for_plot$expr_germline, 0.3))] <- 4
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_germline < quantile(proms_granges_for_plot$expr_germline, 0.3) 
                                   & proms_granges_for_plot$expr_germline > quantile(proms_granges_for_plot$expr_germline, 0.2))] <- 3
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_germline < quantile(proms_granges_for_plot$expr_germline, 0.2) 
                                   & proms_granges_for_plot$expr_germline > quantile(proms_granges_for_plot$expr_germline, 0.1))] <- 2
proms_granges_for_plot$group[which(proms_granges_for_plot$expr_germline < quantile(proms_granges_for_plot$expr_germline, 0.1) 
)] <- 1

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                hypometh_size = proms_granges_for_plot$hypometh_size_human)



quartz(file = "hypometh_size_expr_germline.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = hypometh_size, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), scale = 1.1, from = 0, to = 7500) + 
  labs(x = "hypometh region size", y = "testis expression decile") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size = 8)) + theme(axis.title=element_text(size = 10)) +
  #geom_vline(xintercept = mode_to_plot)
  geom_vline(xintercept = median(proms_df_for_plot$hypometh_size[which(proms_df_for_plot$group1 == 10)]))

dev.off()



####
quartz(file = "hypometh_size_histogram_human.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
hist(tss_with_hypometh$hypometh_size_human, col = alpha("red", 0.6), lty = 0, 
     breaks = 65, freq = FALSE, xlab = "hypomethylated region size (bp)", cex.lab = 1.2, xlim = c(0, 10000),
     yaxt = 'n', xaxt = 'n', main = "", font.main = 1, cex.main = 1.25)
axis(1, at = c(0, 5000, 10000), cex.axis = 1.2)
axis(2, at = c(0, 0.0003), cex.axis = 1.2)
#abline(v = c(40, 80), lty = "longdash", col = rgb(0,0,0,0.7))
dev.off()


























