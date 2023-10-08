getMutationRateDF <- function(mutation_data_frame, mutation_type, region_granges, 
                                   region_type = c("exon", "transcript"), 
                                   use_pseudocount = c(TRUE, FALSE), pseudocount){
  mut_df_restricted <- mutation_data_frame[which(mutation_data_frame$ExonicFunc.refGene %in% mutation_type), ]
  de_novo_mutations <- makeGRangesFromDataFrame(mut_df_restricted, keep.extra.columns = TRUE)
  
  genome(de_novo_mutations) <- "hg19" #based on paper it is hg19
  seqlevelsStyle(de_novo_mutations) <- "ucsc"
  
  de_novo_overlaps <- findOverlaps(region_granges, de_novo_mutations)
  region_granges$de_novo_count <- 0
  mut_by_region <- split(de_novo_mutations[subjectHits(de_novo_overlaps)], queryHits(de_novo_overlaps))
  region_granges$de_novo_count[unique(queryHits(de_novo_overlaps))] <- lengths(mut_by_region)
  
  if (region_type == "exon"){
    regions_by_tx <- split(region_granges, region_granges$tx_id)
    coding_length <- sapply(regions_by_tx, function(xx) sum(width(xx)))
    oe_lof_upper <- constraint[names(regions_by_tx), "oe_lof_upper"]
    de_novo_count <- sapply(regions_by_tx, function(xx) sum(xx$de_novo_count))
    exp_lof_for_filtering <- constraint[names(regions_by_tx), "exp_lof"]
    
    de_novo_count_df <- data.frame(de_novo_count = de_novo_count, cds_length = coding_length, exp_lof = exp_lof_for_filtering,
                                   oe_lof_upper = oe_lof_upper, tx_id = names(regions_by_tx),
                                   stringsAsFactors = FALSE)
    de_novo_count_df$mut_rate <- NA
    de_novo_count_df$mut_rate[which(de_novo_count_df$de_novo_count > 0)] <- de_novo_count_df$de_novo_count[
      which(de_novo_count_df$de_novo_count > 0)]/de_novo_count_df$cds_length[which(de_novo_count_df$de_novo_count > 0)]
    
    if (use_pseudocount == TRUE){
      de_novo_count_df$mut_rate_2 <- (de_novo_count_df$de_novo_count + pseudocount)/de_novo_count_df$cds_length
    }
    regions_for_plot <- de_novo_count_df[which(de_novo_count_df$exp_lof >= 10), ] 
  } 
  else if (region_type == "transcript"){
    region_granges$gene_length <- width(region_granges)
    region_granges$oe_lof_upper <- constraint[names(region_granges), "oe_lof_upper"] #constraint is from the gnomAD paper
    region_granges$exp_lof <- constraint[names(region_granges), "exp_lof"]
    region_granges$mut_rate <- NA
    region_granges$mut_rate[which(region_granges$de_novo_count > 0)] <- region_granges$de_novo_count[
      which(region_granges$de_novo_count > 0)]/region_granges$gene_length[which(region_granges$de_novo_count > 0)]
    region_granges$mut_rate_including_zeros <- region_granges$de_novo_count/region_granges$gene_length
    
    region_granges_for_loeuf <- region_granges[which(region_granges$exp_lof >= 10)]
    
    regions_for_plot <- region_granges_for_loeuf
  }
  
  regions_for_plot$group <- 10
  regions_for_plot$group[which(regions_for_plot$oe_lof_upper < quantile(regions_for_plot$oe_lof_upper, 0.9) 
                               & regions_for_plot$oe_lof_upper > quantile(regions_for_plot$oe_lof_upper, 0.8))] <- 9
  regions_for_plot$group[which(regions_for_plot$oe_lof_upper < quantile(regions_for_plot$oe_lof_upper, 0.8) 
                               & regions_for_plot$oe_lof_upper > quantile(regions_for_plot$oe_lof_upper, 0.7))] <- 8
  regions_for_plot$group[which(regions_for_plot$oe_lof_upper < quantile(regions_for_plot$oe_lof_upper, 0.7) 
                               & regions_for_plot$oe_lof_upper > quantile(regions_for_plot$oe_lof_upper, 0.6))] <- 7
  regions_for_plot$group[which(regions_for_plot$oe_lof_upper < quantile(regions_for_plot$oe_lof_upper, 0.6) 
                               & regions_for_plot$oe_lof_upper > quantile(regions_for_plot$oe_lof_upper, 0.5))] <- 6
  regions_for_plot$group[which(regions_for_plot$oe_lof_upper < quantile(regions_for_plot$oe_lof_upper, 0.5) 
                               & regions_for_plot$oe_lof_upper > quantile(regions_for_plot$oe_lof_upper, 0.4))] <- 5
  regions_for_plot$group[which(regions_for_plot$oe_lof_upper < quantile(regions_for_plot$oe_lof_upper, 0.4) 
                               & regions_for_plot$oe_lof_upper > quantile(regions_for_plot$oe_lof_upper, 0.3))] <- 4
  regions_for_plot$group[which(regions_for_plot$oe_lof_upper < quantile(regions_for_plot$oe_lof_upper, 0.3) 
                               & regions_for_plot$oe_lof_upper > quantile(regions_for_plot$oe_lof_upper, 0.2))] <- 3
  regions_for_plot$group[which(regions_for_plot$oe_lof_upper < quantile(regions_for_plot$oe_lof_upper, 0.2) 
                               & regions_for_plot$oe_lof_upper > quantile(regions_for_plot$oe_lof_upper, 0.1))] <- 2
  regions_for_plot$group[which(regions_for_plot$oe_lof_upper < quantile(regions_for_plot$oe_lof_upper, 0.1) 
  )] <- 1
  regions_for_plot
}

###get objects that will be used for plots
de_novo_df <- read_csv('de_novo_mutations/all_de_novo_mutations_gene4denovo_updated.csv')
all_SNV_types <- c("-", "frameshift substitution", "nonframeshift substitution", 
                   "nonsynonymous SNV", "nonsynonymous SNV;splicing",  "stopgain",  "stoploss",  "synonymous SNV", 
                   "synonymous SNV;splicing") #"-" corresponds to mutations outside coding regions

mut_rate_df_syn_no_introns <- getMutationRateDF(de_novo_df, "synonymous SNV", exons_to_use, 
                                                region_type = "exon", use_pseudocount = FALSE)

mut_rate_df_syn_plus_introns <- getMutationRateDF(de_novo_df, c("synonymous SNV", "-"), tx_all, 
                                                  region_type = "transcript", use_pseudocount = FALSE)
mut_rate_df_all_snv_plus_introns <- getMutationRateDF(de_novo_df, all_SNV_types, tx_all, 
                                                      region_type = "transcript", use_pseudocount = FALSE)



mut_rate_df <- data.frame(group = as.factor(mut_rate_df_syn_no_introns$group), 
                               mut_rate = mut_rate_df_syn_no_introns$mut_rate)

quartz(file = "coding_mut_rate_vs_loeuf.pdf", height = 2.4, width = 1.8, pointsize = 8, type = "pdf")
ggplot(mut_rate_df, aes(x = group, y = mut_rate)) + 
  geom_boxplot(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(0, 0.0025)) + labs(y = "de novo synonymous mutations per bp", x = "LOEUF decile") + 
  theme_classic() + theme(axis.text=element_text(size=7)) + theme(axis.title=element_text(size = 7))
dev.off()

quartz(file = "coding_mut_rate_1vs10_loeuf_deciles.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(density(na.omit(mut_rate_df_syn_no_introns$mut_rate[which(mut_rate_df_syn_no_introns$group == 1)]), 
             from = 0), col = "orange", lwd = 2.5, bty = 'l', xlab = "", main = "", xaxt = 'n', yaxt = 'n')
lines(density(na.omit(mut_rate_df_syn_no_introns$mut_rate[which(mut_rate_df_syn_no_introns$group == 10)]), 
              from = 0), col = rgb(0,0,0,0.7), lwd = 2.5)
axis(1, at = c(0, 0.002, 0.004), labels = c("0", "0.002", "0.004"))
axis(2, at = c(0, 1000))
legend = legend("topright", legend = c("1st LOEUF decile", "10th LOEUF decile"), 
                col = c("orange", rgb(0,0,0,0.7)), lwd = 2.5, bty = 'n')
dev.off()

mod_no_introns <- lm(mut_rate_df_syn_no_introns$mut_rate ~ mut_rate_df_syn_no_introns$oe_lof_upper)
mod_with_introns <- lm(mut_rate_df_all_snv_plus_introns$mut_rate ~ mut_rate_df_all_snv_plus_introns$oe_lof_upper)
mod_with_introns_2 <- lm(mut_rate_df_syn_plus_introns$mut_rate ~ mut_rate_df_syn_plus_introns$oe_lof_upper)
mod_no_introns_without_h3k36me3 <- lm(mut_rate_df_syn_no_introns_non_h3k36me3$mut_rate ~ mut_rate_df_syn_no_introns_non_h3k36me3$oe_lof_upper)

quartz(file = "loeuf_mut_rate_r_squared.pdf", height = 2.5, width = 1.8, pointsize = 8, type = "pdf")
par(mar = c(4, 5.5, 1, 1))
plot(1, 100*summary(mod_no_introns)$adj.r.squared, pch = 19, cex = 1.5, col = alpha("red", 0.62), bty = 'l', 
     xlab = "", ylab = "% explained\nmut rate variance by LOEUF", cex.lab = 1.1, ylim = c(0, 4), xlim = c(0.8, 4.2), xaxt = 'n', yaxt = 'n')
points(2, 100*summary(mod_no_introns_without_h3k36me3)$adj.r.squared, pch = 19, cex = 1.5, col = alpha("red", 0.62))
points(3, 100*summary(mod_with_introns_2)$adj.r.squared, pch = 19, cex = 1.5, col = alpha("red", 0.62))
points(4, 100*summary(mod_with_introns)$adj.r.squared, pch = 19, cex = 1.5, col = alpha("red", 0.62))
axis(2, at = c(0, 2, 4), cex.axis = 1.2)
axis(1, at = c(1,2,3,4), labels = c("synonymous", "synonymous w/out\nH3K36me3 regions", "synonymous\n+ intronic", "all coding SNV\n+intronic"), 
     las = 2, cex.axis = 0.5)
dev.off()




###pseudocount
#now plot the same relationship between de novo mut rate and loeuf across all genes, by using the pseudocount of 0.5 to each gene (so now there's no genes with 0 counts)
mut_rate_df_syn_no_introns <- getMutationRateDF(de_novo_df, "synonymous SNV", exons_to_use, 
                                                region_type = "exon", use_pseudocount = TRUE, 
                                                pseudocount = 1)

mut_rate_df2 <- data.frame(group = as.factor(mut_rate_df_syn_no_introns$group), 
                          mut_rate = mut_rate_df_syn_no_introns$mut_rate_2)

quartz(file = "coding_mut_rate_vs_loeuf_2.pdf", height = 2.4, width = 1.8, pointsize = 8, type = "pdf")
ggplot(mut_rate_df2, aes(x = group, y = mut_rate)) + 
  geom_boxplot(fill = alpha("red", 0.59), color = rgb(0,0,0,0.5), outlier.shape=NA) + 
  coord_cartesian(ylim = c(0, 0.0032)) + labs(y = "de novo synonymous mutations per bp", x = "LOEUF decile") + 
  theme_classic() + theme(axis.text=element_text(size=7)) + theme(axis.title=element_text(size = 7))
dev.off()





