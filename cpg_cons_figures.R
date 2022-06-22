proms_granges_for_plot <- prom_proximal_score_df 
#or for promoter boundary CpGs do this instead: proms_granges_for_plot <- prom_boundary_score_df
proms_granges_for_plot$group <- 4
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.75) 
                                   & proms_granges_for_plot$oe_lof_upper >= quantile(proms_granges_for_plot$oe_lof_upper, 0.5))] <- 3
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.5) 
                                   & proms_granges_for_plot$oe_lof_upper >= quantile(proms_granges_for_plot$oe_lof_upper, 0.25))] <- 2
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper <= quantile(proms_granges_for_plot$oe_lof_upper, 0.25) )] <- 1

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, levels = c(1, 2, 3, 4))




proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                tissue_spec = proms_granges_for_plot$score_cpgs)


quartz(file = "CpGs_vs_loeuf_phyloP_core.pdf", width = 1.8, height = 2.4, pointsize = 8, type = "pdf")
par(mar = c(6, 4, 1, 1)+0.2)
ggplot(proms_df_for_plot, aes(x = group1, y = tissue_spec)) + 
  geom_boxplot(fill = alpha("red", 0.75), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(-1.3, 1.5)) + labs(y = "mean PhyloP \n(100 vertebrates)", x = "LOEUF quartile") + 
  theme(axis.text=element_text(size=10)) + theme_classic()
dev.off()

proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                tissue_spec = proms_granges_for_plot$score_diff)

quartz(file = "CpGs_vs_non_CpGs_phyloP_core.pdf", width = 1.8, height = 2.4, pointsize = 8, type = "pdf")
par(mar = c(6, 4, 1, 1)+0.2)
ggplot(proms_df_for_plot, aes(x = group1, y = tissue_spec)) + 
  geom_boxplot(fill = alpha("red", 0.75), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(-1.15, 0.5)) + labs(y = "mean PhyloP difference\n(CpGs vs non-CpGs)", x = "LOEUF quartile") + 
  theme(axis.text=element_text(size=10)) + theme_classic()
dev.off()

##promoter boundary CpGs now, remember to generate 
##a new "proms_df_for_plot" as described at the beginning of the script
quartz(file = "CpGs_vs_loeuf_phyloP_boundary.pdf", width = 1.8, height = 2.4, pointsize = 8, type = "pdf")
par(mar = c(6, 4, 1, 1)+0.2)
ggplot(proms_df_for_plot, aes(x = group1, y = tissue_spec)) + 
  geom_boxplot(fill = alpha("orange"), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(-1.1, 0.9)) + labs(y = "mean PhyloP \n(100 vertebrates)", x = "LOEUF quartile") + 
  theme(axis.text=element_text(size=10)) + theme_classic()
dev.off()

proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                tissue_spec = proms_granges_for_plot$score_diff)

quartz(file = "CpGs_vs_non_CpGs_phyloP_boundary.pdf", width = 1.8, height = 2.4, pointsize = 8, type = "pdf")
par(mar = c(6, 4, 1, 1)+0.2)
ggplot(proms_df_for_plot, aes(x = group1, y = tissue_spec)) + 
  geom_boxplot(fill = alpha("orange"), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(-1, 0.1)) + labs(y = "mean PhyloP difference\n(CpGs vs non-CpGs)", x = "LOEUF quartile") + 
  theme(axis.text=element_text(size=10)) + theme_classic()
dev.off()









###
proms_granges_for_plot <- proms_good_power_init2
proms_granges_for_plot$group <- NA
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper > quantile(proms_granges_for_plot$oe_lof_upper, 0.9))] <- 3
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.55) 
                                   & proms_granges_for_plot$oe_lof_upper >= quantile(proms_granges_for_plot$oe_lof_upper, 0.45))] <- 2

proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper <= quantile(proms_granges_for_plot$oe_lof_upper, 0.1))] <- 1

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, levels = c(1, 2, 3))

proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                tissue_spec = proms_granges_for_plot$score_cpgs)

proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                tissue_spec = proms_granges_for_plot$score_diff)






####
proms_granges_for_plot <- prom_boundary_score_df
proms_granges_for_plot$group <- 4
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.75) 
                                   & proms_granges_for_plot$oe_lof_upper >= quantile(proms_granges_for_plot$oe_lof_upper, 0.5))] <- 3
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.5) 
                                   & proms_granges_for_plot$oe_lof_upper >= quantile(proms_granges_for_plot$oe_lof_upper, 0.25))] <- 2
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper <= quantile(proms_granges_for_plot$oe_lof_upper, 0.25) )] <- 1

proms_granges_for_plot$group <- factor(proms_granges_for_plot$group, levels = c(1, 2, 3, 4))

proms_df_for_plot <- data.frame(group1 = proms_granges_for_plot$group, 
                                tissue_spec = proms_granges_for_plot$score_cpgs)

quartz(file = "CpGs_vs_loeuf_phyloP_core.pdf", width = 1.8, height = 2.4, pointsize = 8, type = "pdf")
par(mar = c(6, 4, 1, 1)+0.2)
ggplot(proms_df_for_plot, aes(x = group1, y = tissue_spec)) + 
  geom_boxplot(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(-1.1, 0.9)) + labs(y = "mean PhyloP score\n(100 vertebrates)", x = "LOEUF quartile") + 
  theme(axis.text=element_text(size=10)) + theme_classic()
dev.off()



