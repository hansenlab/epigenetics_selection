#####
library(GenomicScores)
gsco <- getGScores("phyloP100way.UCSC.hg19")
#####

###CpGs in promoter boundary (defined as hypomethylated region overlapping TSS excluding the region +/- 500bp)
hypometh_regions_human <- import("~/Downloads/germline_methylation/human_germline_hypometh_regions.bed", 
                                 format = "BED")
genome(hypometh_regions_human) <- "hg19"
hypometh_regions_human <- hypometh_regions_human[which(seqnames(hypometh_regions_human) %in% names(Hsapiens)[1:22])]
seqlevels(hypometh_regions_human) <- seqlevels(hypometh_regions_human)[1:22]

overlaps <- findOverlaps(hypometh_regions_human, tss_with_hypometh) #"tss_with hypometh" granges defined in "cpgs_meth_promoter_boundary.R" script
dup_indices <- unique(queryHits(overlaps)[which(duplicated(queryHits(overlaps)))])
hypometh_regions_human <- hypometh_regions_human[-dup_indices] #this excludes regions that overlap more than one tss

overlaps <- findOverlaps(hypometh_regions_human, tss_with_hypometh)
hypometh_regions_human <- hypometh_regions_human[unique(queryHits(overlaps))]
hypometh_regions_human$tx_id <- tss_with_hypometh$tx_id[subjectHits(overlaps)]
hypometh_regions_human$oe_lof_upper <- tss_with_hypometh$oe_lof_upper[subjectHits(overlaps)]


###
cpgs_with_meth_level <- cpgs_with_meth_level[which(seqnames(cpgs_with_meth_level) %in% names(Hsapiens)[1:22])]
seqlevels(cpgs_with_meth_level) <- seqlevels(cpgs_with_meth_level)[1:22]

overlaps <- findOverlaps(cpgs_with_meth_level, hypometh_regions_human)
cpgs_prom_boundary <- cpgs_with_meth_level[queryHits(overlaps)]
cpgs_prom_boundary$oe_lof_upper <- hypometh_regions_human$oe_lof_upper[subjectHits(overlaps)]
cpgs_prom_boundary$tx_id <- hypometh_regions_human$tx_id[subjectHits(overlaps)]
cpgs_prom_boundary <- cpgs_prom_boundary[-queryHits(findOverlaps(cpgs_prom_boundary, cpgs_proms_good_power))]

#get CpG granges with both the C and the G as separate ranges of width 1
genome(cpgs_prom_boundary) <- "hg19"

cpgs_prom_boundary_2 <- makeGRangesFromDataFrame(data.frame(chr = seqnames(cpgs_prom_boundary), 
                                                       start = start(cpgs_prom_boundary)+1, 
                                                       end = end(cpgs_prom_boundary)+1, 
                                                       meth_level = cpgs_prom_boundary$meth_level, 
                                                       coverage = cpgs_prom_boundary$coverage, 
                                                       tx_id = cpgs_prom_boundary$tx_id), 
                                            keep.extra.columns = TRUE)


cpgs_prom_boundary_both_strands <- unlist(GRangesList(cpgs_prom_boundary, cpgs_prom_boundary_2))
cpg_ranges <- cpgs_prom_boundary_both_strands
cpg_ranges$score <- gscores(gsco, cpg_ranges)$default

prom_boundary_no_cgs <- setdiff(setdiff(hypometh_regions_human, proms_good_power-1500, ignore.strand = TRUE), 
                        cpg_ranges, ignore.strand = TRUE)

overlaps <- findOverlaps(prom_boundary_no_cgs, hypometh_regions_human)
prom_boundary_no_cgs <- prom_boundary_no_cgs[queryHits(overlaps)]
prom_boundary_no_cgs$tx_id <- hypometh_regions_human$tx_id[subjectHits(overlaps)]

prom_boundary_no_cgs$score <- gscores(gsco, prom_boundary_no_cgs)$default

save(prom_boundary_no_cgs, cpg_ranges, file = "promoter_boundary_nucleotides_phyloP.rda")

cpgs_by_tx <- split(cpg_ranges, cpg_ranges$tx_id)
prom_boundary_by_tx <- split(prom_boundary_no_cgs, prom_boundary_no_cgs$tx_id)

cpgs_by_tx <- cpgs_by_tx[names(prom_boundary_by_tx)]
prom_boundary_by_tx <- prom_boundary_by_tx[names(cpgs_by_tx)]

prom_boundary_score_df <- data.frame(tx_id = names(cpgs_by_tx), 
                            score_cpgs = as.numeric(sapply(cpgs_by_tx, function(xx) mean(xx$score))),
                            score_no_cpgs = as.numeric(sapply(prom_boundary_by_tx, function(xx) {
                              gr <- xx
                              widths <- width(gr)
                              weighted_average_numerator <- sum(gr$score*widths)
                              weighted_average_denominator <- sum(widths)
                              weighted_average_numerator/weighted_average_denominator
                            })),
                            stringsAsFactors = FALSE)

prom_boundary_score_df$oe_lof_upper <- sapply(prom_boundary_score_df$tx_id, 
                                              function(xx) proms_good_power$oe_lof_upper[which(proms_good_power$tx_id == xx)])
prom_boundary_score_df$score_diff <- prom_boundary_score_df$score_cpgs - prom_boundary_score_df$score_no_cpgs

###CpGs in proximal promoter (+/- 500bp from TSS)
overlaps <- findOverlaps(cpgs_with_meth_level, proms_good_power-1500) #proms_good_power granges have width 4000 each so we do -1500 to get proximal promoter
cpgs_proms_good_power <- cpgs_with_meth_level[unique(queryHits(overlaps))]
cpgs_proms_good_power$tx_id <- proms_good_power$tx_id[subjectHits(overlaps)]

cpgs_in_proms <- cpgs_proms_good_power #all CpGs that overlap proms (defined as regions +/- 2kb from TSS)
genome(cpgs_in_proms) <- "hg19"

cpgs_in_proms_2 <- makeGRangesFromDataFrame(data.frame(chr = seqnames(cpgs_in_proms), 
                                                       start = start(cpgs_in_proms)+1, 
                                                       end = end(cpgs_in_proms)+1, 
                                                       meth_level = cpgs_in_proms$meth_level, 
                                                       coverage = cpgs_in_proms$coverage, 
                                                       tx_id = cpgs_in_proms$tx_id), 
                                            keep.extra.columns = TRUE)


cpgs_in_proms_both_strands <- unlist(GRangesList(cpgs_in_proms, cpgs_in_proms_2))
cpg_ranges <- cpgs_in_proms_both_strands
cpg_ranges$score <- gscores(gsco, cpg_ranges)$default

prom_proximal_no_cgs <- setdiff(proms_good_power-1500, cpg_ranges, ignore.strand = TRUE)
overlaps <- findOverlaps(prom_proximal_no_cgs, proms_good_power-1500)
prom_proximal_no_cgs$tx_id <- proms_good_power$tx_id[subjectHits(overlaps)] #can do this because all queryHits are unique
prom_proximal_no_cgs$score <- gscores(gsco, prom_proximal_no_cgs)$default


#no_cpg_regions_by_prom <- split(prom_proximal_no_cgs[subjectHits(overlaps)], queryHits(overlaps))
#names(no_cpg_regions_by_prom) <- proms_good_power$tx_id
#no_cpg_regions_by_prom <- endoapply(no_cpg_regions_by_prom, function(xx) {
#  xx$tx_id <- names(xx)
#  xx
#})
save(prom_proximal_no_cgs, cpg_ranges, file = "promoter_proximal_nucleotides_phyloP.rda")

cpgs_by_tx <- split(cpg_ranges, cpg_ranges$tx_id)
prom_proximal_by_tx <- split(prom_proximal_no_cgs, prom_proximal_no_cgs$tx_id)

prom_proximal_by_tx <- prom_proximal_by_tx[names(cpgs_by_tx)]
cpgs_by_tx <- cpgs_by_tx[names(prom_proximal_by_tx)]

prom_proximal_score_df <- data.frame(tx_id = names(cpgs_by_tx), 
                                     score_cpgs = as.numeric(sapply(cpgs_by_tx, function(xx) mean(xx$score))),
                                     score_no_cpgs = as.numeric(sapply(prom_proximal_by_tx, function(xx) {
                                       gr <- xx
                                       widths <- width(gr)
                                       weighted_average_numerator <- sum(gr$score*widths)
                                       weighted_average_denominator <- sum(widths)
                                       weighted_average_numerator/weighted_average_denominator
                                     })),
                                     stringsAsFactors = FALSE)

prom_proximal_score_df$oe_lof_upper <- sapply(prom_proximal_score_df$tx_id, 
                                              function(xx) proms_good_power$oe_lof_upper[which(proms_good_power$tx_id == xx)])
prom_proximal_score_df$score_diff <- prom_proximal_score_df$score_cpgs - prom_proximal_score_df$score_no_cpgs



