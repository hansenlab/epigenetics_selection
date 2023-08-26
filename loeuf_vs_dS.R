###assess relationship between dS and loeuf
library(biomaRt)
###need an older version of biomart because for some reason newest versions don't include dn and ds
ensembl96 <- useEnsembl(biomart = 'ensembl', 
                        dataset = 'hsapiens_gene_ensembl',
                        version = 96)

genes_with_ds_1 <- getBM(attributes = c("ensembl_transcript_id", "external_gene_name",
                             "ptroglodytes_homolog_ensembl_gene", 
                             "ptroglodytes_homolog_ds", 
                             "ptroglodytes_homolog_dn"),
              filters = "ensembl_transcript_id",
              values = de_novo_count_df$tx_id,
              mart = ensembl96)
genes_with_ds <- genes_with_ds_1
genes_with_ds$oe_lof_upper <- constraint[genes_with_ds$ensembl_transcript_id, "oe_lof_upper"]

genes_with_ds_no_duplicates <- genes_with_ds[-which(genes_with_ds$ensembl_transcript_id %in% 
                                     genes_with_ds$ensembl_transcript_id[which(duplicated(genes_with_ds$ensembl_transcript_id))]), ]
genes_with_ds_no_duplicates <- genes_with_ds_no_duplicates[-which(genes_with_ds_no_duplicates$ptroglodytes_homolog_ds == 0), ]


proms_granges_for_plot <- genes_with_ds_no_duplicates
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
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.1, na.rm = TRUE) 
)] <- 1



#plot the relationship between de novo mut rate and loeuf first by excluding genes with 0 de novo synonymous mutations
proms_df_for_plot <- data.frame(group1 = as.factor(proms_granges_for_plot$group), 
                                ds = proms_granges_for_plot$ptroglodytes_homolog_ds)

quartz(file = "loeuf_vs_dS.pdf", height = 2.4, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = group1, y = ds)) + 
  geom_boxplot(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(0, 0.07)) + labs(y = "dS", x = "LOEUF decile") + 
  theme_classic() + theme(axis.text=element_text(size=7)) + theme(axis.title=element_text(size = 7))
dev.off()

quartz(file = "loeuf_vs_dS.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = ds, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), scale = 0.9, from = 0, to = 0.05) + 
  labs(x = "dS (human vs chimp)", y = "LOEUF decile") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size = 8)) + theme(axis.title=element_text(size = 10))
dev.off()


###do the same analysis this time vs mouse
genes_with_ds_2 <- getBM(attributes = c("ensembl_transcript_id", "external_gene_name",
                                      "mmusculus_homolog_ensembl_gene", 
                                      "mmusculus_homolog_ds", 
                                      "mmusculus_homolog_dn"),
                       filters = "ensembl_transcript_id",
                       values = de_novo_count_df$tx_id,
                       mart = ensembl96)
genes_with_ds <- genes_with_ds_2
genes_with_ds$oe_lof_upper <- constraint[genes_with_ds$ensembl_transcript_id, "oe_lof_upper"]

genes_with_ds_no_duplicates <- genes_with_ds[-which(genes_with_ds$ensembl_transcript_id %in% 
                                                      genes_with_ds$ensembl_transcript_id[which(duplicated(genes_with_ds$ensembl_transcript_id))]), ]
genes_with_ds_no_duplicates <- genes_with_ds_no_duplicates[-which(genes_with_ds_no_duplicates$mmusculus_homolog_ds == 0), ]


proms_granges_for_plot <- genes_with_ds_no_duplicates
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
proms_granges_for_plot$group[which(proms_granges_for_plot$oe_lof_upper < quantile(proms_granges_for_plot$oe_lof_upper, 0.1, na.rm = TRUE) 
)] <- 1



#plot the relationship between de novo mut rate and loeuf first by excluding genes with 0 de novo synonymous mutations
proms_df_for_plot <- data.frame(group1 = as.factor(proms_granges_for_plot$group), 
                                ds = proms_granges_for_plot$mmusculus_homolog_ds)



quartz(file = "loeuf_vs_dS_2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
values_high_decile <- proms_df_for_plot$ds[which(proms_df_for_plot$group1 == 1)]
values_high_decile <- values_high_decile[-which(is.na(values_high_decile))]
mode_to_plot <- density(values_high_decile)$x[which.max(density(values_high_decile)$y)]
ggplot(proms_df_for_plot, aes(x = ds, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), scale = 1.2, from = 0, to = 1.2) + 
  labs(x = "dS (human vs mouse)", y = "LOEUF decile") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size = 8)) + theme(axis.title=element_text(size = 10)) +  geom_vline(xintercept = mode_to_plot) 
dev.off()
