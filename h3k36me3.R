###
#germline_h3K36me3 <- read_tsv('h3k36me3_sperm_peaks.narrowPeak', col_names = FALSE) #H3K36me3 is a broad mark so I use the broad peaks instead
germline_h3K36me3 <- read_tsv('h3k36me3_sperm_peaks.broadPeak', col_names = FALSE)
colnames(germline_h3K36me3)[c(1,2,3,7)] <- c("chr", "start", "end","signalValue")
h3k36me3 <- makeGRangesFromDataFrame(germline_h3K36me3, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)

###create data frame used for visualization. 
###Needs the "exons_to_use" object generated in the "lof_intolerance_de_novo_variants.R" script 
overlaps <- findOverlaps(exons_to_use, h3k36me3)

k36me3_df <- data.frame(tx_id = names(exons_by_tx), 
                                has_h3k36me3_peak = rep("no", length(names(exons_by_tx))),  
                                stringsAsFactors = FALSE)
rownames(k36me3_df) <- k36me3_df$tx_id
k36me3_df[unique(exons_to_use$tx_id[queryHits(overlaps)]), "has_h3k36me3_peak"] <- "yes"
k36me3_df$loeuf <- constraint[rownames(k36me3_df), "oe_lof_upper"]
k36me3_df$exp_lof <- constraint[rownames(k36me3_df), "exp_lof"]
k36me3_df <- k36me3_df[order(k36me3_df$loeuf), ]
k36me3_df <- k36me3_df[which(k36me3_df$exp_lof >= 10), ]


###
h3k36me3_vs_loeuf_vec <- c(getPercentOfGenesWithCodingPeak(k36me3_df$has_h3k36me3_peak, k36me3_df$loeuf, 0.1, "less"), 
                       getPercentOfGenesWithCodingPeak(k36me3_df$has_h3k36me3_peak, k36me3_df$loeuf, c(0.1, 0.2), "in between"),
                       getPercentOfGenesWithCodingPeak(k36me3_df$has_h3k36me3_peak, k36me3_df$loeuf, c(0.2, 0.3), "in between"),
                       getPercentOfGenesWithCodingPeak(k36me3_df$has_h3k36me3_peak, k36me3_df$loeuf, c(0.3, 0.4), "in between"),
                       getPercentOfGenesWithCodingPeak(k36me3_df$has_h3k36me3_peak, k36me3_df$loeuf, c(0.4, 0.5), "in between"),
                       getPercentOfGenesWithCodingPeak(k36me3_df$has_h3k36me3_peak, k36me3_df$loeuf, c(0.5, 0.6), "in between"),
                       getPercentOfGenesWithCodingPeak(k36me3_df$has_h3k36me3_peak, k36me3_df$loeuf, c(0.6, 0.7), "in between"),
                       getPercentOfGenesWithCodingPeak(k36me3_df$has_h3k36me3_peak, k36me3_df$loeuf, c(0.7, 0.8), "in between"),
                       getPercentOfGenesWithCodingPeak(k36me3_df$has_h3k36me3_peak, k36me3_df$loeuf, c(0.8, 0.9), "in between"),
                       getPercentOfGenesWithCodingPeak(k36me3_df$has_h3k36me3_peak, k36me3_df$loeuf, 0.9, "greater"))

quartz(file = "h3k36me3_vs_loeuf_coding_seq.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(h3k36me3_vs_loeuf_vec, pch = 19, cex = 1.2, col = alpha("red", 0.75), bty = 'l', ylim = c(0.22, 0.6),
     main = "coding sequence", font.main = 1, ylab = "% genes w/ H3K36me3 peak", xlab = "LOEUF decile", 
     xaxt = 'n', yaxt = 'n')
axis(1, at = 1:10, cex.axis = 0.7)
axis(2, at = c(0.25, 0.55))
dev.off()




