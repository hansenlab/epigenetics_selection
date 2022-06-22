de_novo_df <- read_tsv('de_novo_mutations/de_novo_mutations_gene4denovo.txt')
de_novo_df <- de_novo_df[which(de_novo_df$ExonicFunc.refGene %in% c("synonymous SNV")), ]
#de_novo_df <- de_novo_df[which(de_novo_df$ExonicFunc.refGene == "synonymous SNV" | 
#                                 (de_novo_df$ExonicFunc.refGene == "nonsynonymous SNV" & de_novo_df$Phenotype == "Control")), ]
de_novo_df <- de_novo_df[, c("Chr", "Start", "End", "Ref", "Alt", "ExonicFunc.refGene", "Phenotype", "Study")]
de_novo_mutations <- makeGRangesFromDataFrame(de_novo_df, keep.extra.columns = TRUE)
genome(de_novo_mutations) <- "hg19" #based on paper it is hg19
seqlevelsStyle(de_novo_mutations) <- "ucsc"

###do the analysis it using all canonical coding sequences
constraint <- read_tsv('constraint.txt') #the "constraint.txt" file is from gnomAD 
constraint <- as(constraint, "DataFrame")
rownames(constraint) <- constraint$transcript

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
all_coding_exons <- cdsBy(edb, filter = SeqNameFilter(c(1:22)), columns = "tx_id")
all_coding_exons <- unlist(all_coding_exons)
genome(seqinfo(all_coding_exons)) <- "hg19"
seqlevelsStyle(all_coding_exons) <- "ucsc"
exons_to_use <- all_coding_exons[which(all_coding_exons$tx_id %in% constraint$transcript[which(constraint$canonical == TRUE)])]

#
de_novo_overlaps <- findOverlaps(exons_to_use, de_novo_mutations)
exons_to_use$de_novo_count <- 0
mut_by_exon <- split(de_novo_mutations[subjectHits(de_novo_overlaps)], queryHits(de_novo_overlaps))
exons_to_use$de_novo_count[unique(queryHits(de_novo_overlaps))] <- lengths(mut_by_exon)

exons_by_tx <- split(exons_to_use, exons_to_use$tx_id)

coding_length <- sapply(exons_by_tx, function(xx) sum(width(xx)))
oe_lof_upper <- constraint[names(exons_by_tx), "oe_lof_upper"]
de_novo_count <- sapply(exons_by_tx, function(xx) sum(xx$de_novo_count))
exp_lof_for_filtering <- constraint[names(exons_by_tx), "exp_lof"]

de_novo_count_df <- data.frame(de_novo_count = de_novo_count, cds_length = coding_length, exp_lof = exp_lof_for_filtering,
                               oe_lof_upper = oe_lof_upper, tx_id = names(exons_by_tx),
                               stringsAsFactors = FALSE)
de_novo_count_df$mut_rate <- NA
de_novo_count_df$mut_rate[which(de_novo_count_df$de_novo_count > 0)] <- de_novo_count_df$de_novo_count[
  which(de_novo_count_df$de_novo_count > 0)]/de_novo_count_df$cds_length[which(de_novo_count_df$de_novo_count > 0)]

de_novo_count_df$mut_rate_2 <- (de_novo_count_df$de_novo_count + 1)/de_novo_count_df$cds_length

#de_novo_count_df <- de_novo_count_df[-which(is.na(de_novo_count_df$oe_lof_upper)), ]

###For plotting and doing the formal statistics I should restrict to trascripts with exp_lof >= 10
###since otherwise loeuf estimates are not reliable
proms_granges_for_plot <- de_novo_count_df[which(de_novo_count_df$exp_lof >= 10), ] 
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



#plot the relationship between de novo mut rate and loeuf first by excluding genes with 0 de novo synonymous mutations
proms_df_for_plot <- data.frame(group1 = as.factor(proms_granges_for_plot$group), 
                                mut_rate = proms_granges_for_plot$mut_rate)

quartz(file = "coding_mut_rate_vs_loeuf.pdf", height = 2.4, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = group1, y = mut_rate)) + 
  geom_boxplot(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(0, 0.0025)) + labs(y = "de novo synonymous mutations per bp", x = "LOEUF decile") + 
  theme_classic() + theme(axis.text=element_text(size=7)) + theme(axis.title=element_text(size = 7))
dev.off()

#now plot the same relationship between de novo mut rate and loeuf across all genes, by using the pseudocount of 0.5 to each gene (so now there's no genes with 0 counts)
proms_df_for_plot <- data.frame(group1 = as.factor(proms_granges_for_plot$group), 
                                mut_rate = proms_granges_for_plot$mut_rate_2)

quartz(file = "coding_mut_rate_vs_loeuf_2.pdf", height = 2.4, width = 1.8, pointsize = 8, type = "pdf")
ggplot(proms_df_for_plot, aes(x = group1, y = mut_rate)) + 
  geom_boxplot(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(0, 0.0025)) + labs(y = "de novo synonymous mutations per bp", x = "LOEUF decile") + 
  theme_classic() + theme(axis.text=element_text(size=7)) + theme(axis.title=element_text(size = 7))
dev.off()

###ridgeplot instead
proms_df_for_plot <- data.frame(group1 = as.factor(proms_granges_for_plot$group), 
                                mut_rate = proms_granges_for_plot$mut_rate)

quartz(file = "coding_mut_rate_vs_loeuf_ridgeplot.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
values_high_decile <- proms_df_for_plot$mut_rate[which(proms_df_for_plot$group1 == 10)]
values_high_decile <- values_high_decile[-which(is.na(values_high_decile))]
mode_to_plot <- density(values_high_decile)$x[which.max(density(values_high_decile)$y)]
ggplot(proms_df_for_plot, aes(x = mut_rate, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), scale = 1, from = 0, to = 0.003) + 
  labs(x = "de novo synonymous mutations per bp", y = "LOEUF decile") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size = 8)) + theme(axis.title=element_text(size = 10)) + 
  geom_vline(xintercept = mode_to_plot)
dev.off()

proms_df_for_plot <- data.frame(group1 = as.factor(proms_granges_for_plot$group), 
                                mut_rate = proms_granges_for_plot$mut_rate_2)

quartz(file = "coding_mut_rate_vs_loeuf_ridgeplot_2.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
values_high_decile <- proms_df_for_plot$mut_rate[which(proms_df_for_plot$group1 == 10)]
values_high_decile <- values_high_decile[-which(is.na(values_high_decile))]
mode_to_plot <- density(values_high_decile)$x[which.max(density(values_high_decile)$y)]
ggplot(proms_df_for_plot, aes(x = mut_rate, y = group1)) + 
  geom_density_ridges2(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), scale = 1, from = 0, to = 0.003) + 
  labs(x = "de novo synonymous mutations per bp", y = "LOEUF decile") + 
  theme_ridges(center = TRUE, grid = FALSE) + 
  theme(axis.text=element_text(size = 8)) + theme(axis.title=element_text(size = 10)) + 
  geom_vline(xintercept = mode_to_plot)
dev.off()











  
  
  





